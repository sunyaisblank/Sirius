"""Diagnose original retained controls or their exact captured Metal source.

This deliberately bounded producer never writes a qualification receipt. The
runner revision and the reused executable's source revision remain distinct.
"""

from bisect import bisect_right
import hashlib
from pathlib import Path
import re
import struct
import argparse
import json
import os
import resource
import signal
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
import zipfile
from contextlib import ExitStack
import platform

from runtime_identity import provider_inputs, selected_inventory
from macos_metal_probe import FORCED_INLINE_COUNT, MSL_BYTES, MSL_SHA, unforced_msl


def macho_modules(path):
    """Return (metadata, payloads), keyed by qualified retained shader name.

    Accept only a thin little-endian arm64 executable containing the frozen
    33 retained-program and two retained-owner modules. Bounds come from its
    sections and defined symbols, never a scan for SPIR-V magic. This checks
    framing and instruction bounds, not full SPIR-V semantic validity.
    """
    data = Path(path).read_bytes()

    def require(condition, message):
        if not condition:
            raise ValueError("retained Mach-O: " + message)

    def unpack(fmt, offset, limit=None):
        size = struct.calcsize(fmt)
        limit = len(data) if limit is None else limit
        require(0 <= offset <= limit <= len(data) and size <= limit - offset,
                "truncated or out-of-bounds structure")
        return struct.unpack_from(fmt, data, offset)

    header = unpack("<8I", 0)
    require(header[0] == 0xFEEDFACF and header[1] == 0x0100000C and header[3] == 2,
            "expected a thin little-endian arm64 executable")
    commands_end = 32 + header[5]
    require(0 < header[4] <= header[5] // 8 and commands_end <= len(data),
            "invalid load-command bounds")
    sections = []
    symtab = None
    offset = 32
    for _ in range(header[4]):
        command, size = unpack("<II", offset, commands_end)
        require(size >= 8 and size % 8 == 0 and size <= commands_end - offset,
                "invalid load-command size")
        if command == 0x19:  # LC_SEGMENT_64
            segment = unpack("<II16sQQQQiiII", offset, offset + size)
            _, _, segment_name, address, vm_size, file_offset, file_size, _, protection, count, _ = segment
            require(size == 72 + count * 80 and address + vm_size <= 1 << 64,
                    "invalid segment or section count")
            require(file_offset <= len(data) and file_size <= len(data) - file_offset,
                    "segment extends beyond the file")
            for index in range(count):
                section = unpack("<16s16sQQ8I", offset + 72 + index * 80, offset + size)
                name, owner, start, length, position, _, _, _, flags, _, _, _ = section
                require(owner == segment_name and address <= start <= address + vm_size
                        and length <= address + vm_size - start,
                        "section lies outside its segment")
                file_backed = flags & 0xFF not in {1, 0xC, 0x12}
                if file_backed:
                    require(file_offset <= position <= file_offset + file_size
                            and length <= file_offset + file_size - position,
                            "section lies outside its file-backed segment")
                sections.append({"name": name.rstrip(b"\0"),
                                 "segment": owner.rstrip(b"\0"), "address": start,
                                 "bytes": length, "file_offset": position,
                                 "protection": protection, "file_backed": file_backed,
                                 "type": flags & 0xFF})
        elif command == 2:  # LC_SYMTAB
            require(symtab is None and size == 24, "missing or duplicate symbol-table command")
            symtab = unpack("<6I", offset, offset + size)[2:]
        offset += size
    require(offset == commands_end and symtab is not None, "incomplete load commands or symbol table")
    symoff, symbol_count, stroff, string_bytes = symtab
    require(0 < symbol_count <= (len(data) - symoff) // 16 and 0 <= symoff <= len(data),
            "symbol table extends beyond the file")
    require(0 < string_bytes <= len(data) - stroff and 0 <= stroff < len(data),
            "string table extends beyond the file")
    string_end = stroff + string_bytes
    symbols = []
    addresses = {}
    for index in range(symbol_count):
        strx, kind, section_number, _, value = unpack("<IBBHQ", symoff + index * 16)
        if kind & 0xE0 or kind & 0x0E != 0x0E:  # Debug or non-N_SECT symbol.
            continue
        require(1 <= section_number <= len(sections) and strx < string_bytes,
                "invalid defined symbol section or string offset")
        end = data.find(b"\0", stroff + strx, string_end)
        require(end >= 0, "unterminated symbol name")
        name = data[stroff + strx:end]
        symbols.append((name, section_number, value))
        addresses.setdefault(section_number, set()).add(value)
    addresses = {number: sorted(values) for number, values in addresses.items()}
    pattern = re.compile(rb"__ZN6sirius7backend(16retained_program|14retained_owner)(\d+)(k\w+Shader)E")
    expected = {"retained_program::k" + stage + mode + "Shader"
                for stage in ("Camera", "Transport", "Endpoint", "Dense", "Initialize", "RayCamera", "DopriPhase")
                for mode in ("", "Fp64", "Portable", "PortableFp64")}
    expected.update("retained_program::k" + stage + "FmaShader" for stage in ("Transport", "Endpoint"))
    expected.update("retained_program::k" + stage + "PortableNormalSumShader"
                    for stage in ("Transport", "Endpoint", "DopriPhase"))
    expected.update("retained_owner::k" + stage + "Shader" for stage in ("Transport", "Endpoint"))
    metadata, payloads = {}, {}
    for symbol, number, address in symbols:
        match = pattern.fullmatch(symbol)
        if not match:
            continue
        scope, encoded_length, variable = match.groups()
        require(len(variable) == int(encoded_length), "invalid shader symbol length")
        name = scope[2:].decode("ascii") + "::" + variable.decode("ascii")
        require(name in expected and name not in metadata, "unexpected or duplicate shader symbol")
        section = sections[number - 1]
        require(section["segment"] == b"__TEXT" and section["name"] == b"__const"
                and section["file_backed"] and section["type"] == 0
                and section["protection"] & 1 and not section["protection"] & 2,
                "shader is not in allocated read-only constant storage")
        section_end = section["address"] + section["bytes"]
        require(section["address"] <= address < section_end and address % 4 == 0,
                "shader symbol is outside or misaligned in its section")
        following = addresses[number]
        next_index = bisect_right(following, address)
        end = min(section_end, following[next_index] if next_index < len(following) else section_end)
        start = section["file_offset"] + address - section["address"]
        bound = start + end - address
        magic, version, _, id_bound, schema = unpack("<5I", start, bound)
        require(magic == 0x07230203 and version == 0x00010500 and id_bound > 0 and schema == 0,
                "invalid frozen SPIR-V header")
        cursor, instructions, last_opcode = start + 20, 0, None
        while cursor + 4 <= bound:
            word, = unpack("<I", cursor, bound)
            words, opcode = word >> 16, word & 0xFFFF
            if words == 0:
                break
            require(words * 4 <= bound - cursor, "SPIR-V instruction exceeds its symbol boundary")
            cursor += words * 4
            instructions += 1
            last_opcode = opcode
        require(last_opcode == 56 and bound - cursor <= 15 and not any(data[cursor:bound]),
                "SPIR-V has no final function end or has invalid trailing padding")
        payload = data[start:cursor]
        metadata[name] = {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest(),
                          "file_offset": start, "virtual_address": address,
                          "symbol": symbol.decode("ascii"), "section": "__TEXT,__const",
                          "padding_bytes": bound - cursor, "instruction_count": instructions,
                          "spirv_version": version, "id_bound": id_bound}
        payloads[name] = payload
    require(metadata.keys() == expected, "frozen 35-module symbol inventory differs")
    return metadata, payloads


ROOT = Path(__file__).resolve().parent.parent
SOURCE = "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc"
ARCHIVE_SHA = "973d6ada66511b80d2dcf76e0c95a3b38367c32e08a0f16a62c20e0c0a90974d"
BACKEND_SHA = "682b2ef3c91fe2ea46ab97a18096e509cab5d5eb3805a5b886e63acbe00577ad"
APP_SHA = "af691acddbe15b0a3a9a18113b771bd783d280713cc4176f68d1693e3e9aebe1"
GATE_SHA = "68172518e886aa428bbd88420eeed3a05d7ab33e046eea8ff48000f52a66fc8f"
CASES = (
    "PhysicalInitializationRetainsTheHamiltonianResidual",
    "ProjectedEndpointsKeepPhysicalColumnsAndRetainedContinuation",
    "JointRkStagesRetainCriticalIncrementsAndEmbeddedError",
    "DenseSegmentsPreserveSmallCovariantArrivalDerivatives",
)
LIMITS = {"seconds": 300, "rss_kib": 4 * 1024**2, "reserve_kib": 2 * 1024**2}
STOP = None


def require(condition, message):
    if not condition:
        raise ValueError(message)


def write(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def identity(path):
    path = Path(path)
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    after = path.stat()
    require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns),
            f"input changed while hashing: {path}")
    return {"bytes": after.st_size, "sha256": digest.hexdigest()}


def source_state():
    command = lambda *args: subprocess.check_output(["git", *args], cwd=ROOT, text=True)
    revision = command("rev-parse", "HEAD").strip()
    status = command("status", "--porcelain")
    require(not status, "diagnostic requires a clean committed runner checkout")
    changes = set(command("diff", "--name-only", SOURCE, revision, "--").splitlines())
    require(changes <= {".github/workflows/ci.yml", "scripts/diagnose-macos-retained.py",
                        "scripts/governance/build_policy.py", "scripts/macos_metal_probe.py",
                        "scripts/macos_moltenvk_integration.py", "scripts/macos_vulkan_identity.py",
                        "scripts/validate-native-runtime.py", "scripts/macos_moltenvk_provider.py",
                        "README.md"},
            "reused executable's product/scientific source differs from the runner checkout")
    return {"runner_revision": revision, "candidate_revision": SOURCE, "status": status}


def prepare(archive, output):
    require(identity(archive) == {"bytes": 57718924, "sha256": ARCHIVE_SHA},
            "native export differs from pinned artifact 11778215779")
    stage = ROOT / "bin/macos/retained-diagnostic"
    require(not stage.exists(), "diagnostic product directory already exists")
    stage.mkdir(parents=True)
    names = ["native_build_gate.json", "macos-build.json", "alignment_receipt.json",
             "native-build-tested-sirius_backend_tests", "qualification-sirius.bin"]
    with zipfile.ZipFile(archive) as bundle:
        require(len(bundle.infolist()) == 40 and len(set(bundle.namelist())) == 40,
                "native export inventory differs")
        for name in names:
            (stage / name).write_bytes(bundle.read(name))
    app, backend = stage / names[-1], stage / names[-2]
    for path, sha in [(app, APP_SHA), (backend, BACKEND_SHA), (stage / names[0], GATE_SHA)]:
        require(identity(path)["sha256"] == sha, f"pinned payload differs: {path.name}")
    gate = json.loads((stage / names[0]).read_text())
    require(gate["source"] == {"clean": True, "revision": SOURCE}
            and gate["status"] == "passed" and gate["kind"] == "sirius-native-build-gate",
            "native export lacks its original exact-source build gate")
    for path, pin in [(backend, gate["tested_artifacts"]["sirius_backend_tests"]),
                      (app, gate["product_artifacts"]["sirius"]),
                      (stage / names[2], gate["product_artifacts"]["alignment_receipt"])]:
        require(identity(path) == {key: pin[key] for key in ["bytes", "sha256"]}, path.name)
    require(json.loads((stage / names[1]).read_text())["source_revision"] == SOURCE,
            "build record source differs")
    app.chmod(0o755)
    backend.chmod(0o755)
    modules, payloads = macho_modules(backend)
    (stage / "embedded").mkdir()
    for name, payload in payloads.items():
        (stage / "embedded" / (name.replace("::", "-") + ".spv")).write_bytes(payload)
    write(output / "embedded-modules.json", modules)
    write(output / "preparation.json", {
        "source": source_state(), "archive": {"artifact_id": 11778215779, **identity(archive)},
        "executables": {"app": identity(app), "backend": identity(backend)},
        "embedded_named_arrays": len(modules), "embedded_distinct_payloads": len(set(p["sha256"] for p in modules.values())),
        "limits": LIMITS, "qualification_claimed": False,
        "scope": "Exact existing macOS build export, unchanged original fixtures and embedded arrays; no rebuild or runtime admission."})
    return stage, app, backend, modules


def parse_ps(text):
    rows = {}
    for line in text.splitlines():
        fields = line.split(None, 5)
        require(len(fields) == 6, "malformed ps row")
        pid, parent, group, rss = map(int, fields[:4])
        require(pid >= 0 and parent >= 0 and group >= 0 and rss >= 0 and fields[4] and fields[5],
                "invalid ps process fields")
        require(pid not in rows, "duplicate ps process")
        rows[pid] = {"pid": pid, "ppid": parent, "pgid": group, "rss_kib": rss,
                     "state": fields[4], "lstart": fields[5]}
    require(rows, "empty process census")
    return rows


def memory_sample(text, sysctl_text):
    """A sampled Darwin estimate, not a proven Linux MemAvailable equivalent.

    vm_stat's printed free excludes speculative pages. Add only the disjoint
    inactive queue; do not add active, wired, compressed or overlapping pools.
    Inactive anonymous pages may need compression, so retain raw counters and
    require normal kernel pressure as well as the unchanged 2-GiB reserve.
    References: apple-oss-distributions/xnu, xnu-11417.140.69/vm_page.h;
    giampaolo/psutil, 01bb959/psutil/arch/osx/mem.c (raw free + inactive).
    """
    pages = re.findall(r"page size of (\d+) bytes", text)
    require(len(pages) == 1, "unreadable vm_stat page size")
    page_size = int(pages[0])
    require(page_size >= 1024 and page_size & (page_size - 1) == 0, "invalid VM page size")
    counters = {}
    for line in text.splitlines()[1:]:
        match = re.fullmatch(r"([^:]+):\s+(\d+)\.", line.strip())
        require(match is not None and match[1] not in counters, "malformed or duplicate vm_stat counter")
        counters[match[1]] = int(match[2])
    required = ["Pages free", "Pages inactive", "Swapins", "Swapouts", "Compressions",
                "Decompressions", "Pages occupied by compressor", "Pages stored in compressor"]
    require(all(name in counters for name in required), "required VM counters missing")
    values = sysctl_text.splitlines()
    require(len(values) == 2 and all(re.fullmatch(r"\d+", value) for value in values), "unreadable pressure/physical-memory query")
    pressure, physical_bytes = map(int, values)
    available = (counters["Pages free"] + counters["Pages inactive"]) * page_size
    require(physical_bytes > 0 and available <= physical_bytes, "VM estimate exceeds physical memory")
    return {"available_kib": available // 1024, "free_kib": counters["Pages free"] * page_size // 1024,
            "inactive_kib": counters["Pages inactive"] * page_size // 1024,
            "pressure_level": pressure, "physical_bytes": physical_bytes,
            "page_size": page_size, "counters": counters, "raw_vm_stat": text}


def process_census():
    environment = dict(os.environ, LC_ALL="C")
    processes = subprocess.check_output(
        ["/bin/ps", "-axo", "pid=,ppid=,pgid=,rss=,state=,lstart="],
        text=True, env=environment, timeout=5)
    return parse_ps(processes)


def memory_self_test():
    # Regression for the observed 129-MiB-free prelaunch refusal. Only inactive
    # pages supplement free; enormous other counters must not inflate admission.
    fields = {"Pages free": 8256, "Pages inactive": 122816, "Swapins": 0,
              "Swapouts": 0, "Compressions": 900000, "Decompressions": 700000,
              "Pages occupied by compressor": 200000, "Pages stored in compressor": 400000,
              "Pages active": 900000, "Pages speculative": 200000}
    def fixture(values):
        return ("Mach Virtual Memory Statistics: (page size of 16384 bytes)\n"
                + "".join(f"{name}: {value}.\n" for name, value in values.items()))
    text = fixture(fields)
    sample = memory_sample(text, "1\n17179869184\n")
    require(sample["free_kib"] == 132096 and sample["available_kib"] == 2097152,
            "available-memory boundary self-test failed")
    require(memory_sample(fixture({**fields, "Pages inactive": 122815}), "1\n17179869184\n")["available_kib"] < LIMITS["reserve_kib"],
            "below-reserve sample admitted")
    for pressure in [0, 2, 4, 99]:
        require(memory_sample(text, f"{pressure}\n17179869184\n")["pressure_level"] != 1,
                "non-normal pressure normalized")
    for raw, query in [(text + "Pages inactive: 1.\n", "1\n17179869184\n"),
                       (text.replace("Pages free: 8256.", "Pages free: -1."), "1\n17179869184\n"),
                       (text.replace("Swapouts: 0.\n", ""), "1\n17179869184\n"),
                       (text, "1\n"), (text, "1\n1024\n")]:
        try:
            memory_sample(raw, query)
        except ValueError:
            continue
        raise ValueError("malformed memory observation accepted")


def snapshot():
    environment = dict(os.environ, LC_ALL="C")
    processes = process_census()
    memory = subprocess.check_output(["/usr/bin/vm_stat"], text=True, env=environment, timeout=5)
    pressure = subprocess.check_output(
        ["/usr/sbin/sysctl", "-n", "kern.memorystatus_vm_pressure_level", "hw.memsize"],
        text=True, env=environment, timeout=5)
    return processes, memory_sample(memory, pressure)


def child_limits():
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024**2, 64 * 1024**2))


def run_owned(argv, environment, folder, inputs, *, seconds=300, source_snapshot=None):
    """Keep an unreaped child as the killpg anchor; never signal saved bare PIDs.

    ps lstart is a conservative observation token, not a unique kernel birth
    handle. Escaped or unresolved observed children make containment incomplete.
    """
    require(sys.platform == "darwin", "Darwin process observation required")
    require(signal.getsignal(signal.SIGCHLD) == signal.SIG_DFL, "SIGCHLD must not auto-reap children")
    require(STOP is None, "diagnostic interrupted before the next owner launch")
    folder.mkdir()
    before = {str(path): identity(path) for path in inputs}
    write(folder / "inputs-before.json", before)
    observe_source = source_state if source_snapshot is None else source_snapshot
    source = observe_source()
    owner = {"status": "running", "argv": argv, "source": source,
             "limits": {**LIMITS, "seconds": seconds}, "qualification_claimed": False,
             "reserve_scope": "Darwin sampled available estimate: printed free + inactive; kernel pressure must be normal (1). This is not a proven no-swap bound. RSS includes controller and observed group/tree, not global Metal compiler services.",
             "birth_scope": "Observed ps PID/lstart tokens; second-resolution ambiguity remains incomplete.",
             "sampling_interval_seconds": 0.1, "observed_births": [], "remaining_observed": [],
             "signal_observations": [],
             "environment": {key: environment.get(key) for key in ["VK_DRIVER_FILES", "DYLD_LIBRARY_PATH", "SIRIUS_VULKAN_DEVICE", "MVK_CONFIG_SHADER_DUMP_DIR", "MVK_CONFIG_DEBUG", "MVK_CONFIG_TRACE_VULKAN_CALLS"]}}
    child, reaped, reason, code, seen, escaped = None, False, None, None, {}, False
    group_drained = False
    retirement_census_complete = False
    peak, times, cleanup_errors = 0, [], []
    started = time.monotonic()
    def interrupted(number, _):
        global STOP
        STOP = number
    handlers = {}
    def observe(rows):
        nonlocal escaped
        members = {pid for pid, row in rows.items() if row["pgid"] == child.pid}
        members.add(child.pid)
        changed = True
        while changed:
            more = {pid for pid, row in rows.items() if row["ppid"] in members}
            changed = not more <= members
            members |= more
        selected = [rows[pid] for pid in members if pid in rows]
        for row in selected:
            seen[(row["pid"], row["lstart"])] = row
            escaped |= row["pgid"] != child.pid
        return selected
    def drained(rows, members):
        return (child.pid in rows and rows[child.pid]["pgid"] == child.pid
                and rows[child.pid]["state"].startswith("Z")
                and not any(not row["state"].startswith("Z") for row in members))
    def group_signal(number):
        require(child is not None and not reaped, "cannot signal a reaped or absent group anchor")
        try:
            os.killpg(child.pid, number)
        except ProcessLookupError:
            pass
        except PermissionError:
            # Darwin can return EPERM when the last live member exits before
            # killpg. Accept that race only after observing the unreaped zombie
            # anchor and no live group/tree members. Unknown or live targets
            # still fail containment; never retry against a saved bare PID.
            rows = process_census()
            require(drained(rows, observe(rows)), "group signal denied with live or unresolved targets")
            owner["signal_observations"].append({"signal": number, "result": "EPERM; anchored group observed drained"})
    with ExitStack() as files:
        try:
            for number in [signal.SIGINT, signal.SIGTERM, signal.SIGHUP]:
                handlers[number] = signal.signal(number, interrupted)
            stdout = files.enter_context((folder / "stdout.log").open("wb"))
            stderr = files.enter_context((folder / "stderr.log").open("wb"))
            samples = files.enter_context((folder / "samples.jsonl").open("w"))
            write(folder / "owner.json", owner)
            rows, memory = snapshot()
            owner["prelaunch_memory"] = memory
            owner["controller"] = rows[os.getpid()]
            require(memory["available_kib"] >= LIMITS["reserve_kib"], "prelaunch available-memory reserve below 2 GiB")
            require(memory["pressure_level"] == 1, "prelaunch kernel memory pressure is not normal")
            child = subprocess.Popen(argv, cwd=folder, env=environment, stdout=stdout, stderr=stderr,
                                     start_new_session=True, preexec_fn=child_limits)
            owner["pid"] = child.pid
            write(folder / "owner.json", owner)
            while True:
                rows, memory = snapshot()
                members = observe(rows)
                elapsed = time.monotonic() - started
                rss = sum(row["rss_kib"] for row in members) + rows[os.getpid()]["rss_kib"]
                peak = max(peak, rss)
                times.append(elapsed)
                samples.write(json.dumps({"elapsed_seconds": elapsed, "rss_kib": rss,
                                          "memory": memory, "members": members}) + "\n")
                samples.flush()
                if STOP is not None: reason = f"signal {STOP}"
                elif elapsed >= seconds: reason = "deadline"
                elif rss > LIMITS["rss_kib"]: reason = "RSS limit"
                elif memory["available_kib"] < LIMITS["reserve_kib"]: reason = "available-memory reserve"
                elif memory["pressure_level"] != 1: reason = "kernel memory pressure is not normal"
                elif escaped: reason = "observed descendant escaped launch group"
                elif child.pid not in rows: reason = "unreaped leader missing from observation"
                if reason is not None or rows[child.pid]["state"].startswith("Z"):
                    break
                time.sleep(0.1)
        except Exception as error:
            reason = f"observation incomplete: {type(error).__name__}: {error}"
        finally:
            if child is not None:
                # Popen.poll/send_signal/terminate/kill/communicate would reap the
                # anchor. Only killpg is used until all group signaling is over.
                try:
                    rows = process_census()
                    members = observe(rows)
                    if not drained(rows, members):
                        group_signal(signal.SIGTERM)
                    deadline = time.monotonic() + 1
                    while True:
                        rows = process_census()
                        members = observe(rows)
                        group_drained = drained(rows, members)
                        if group_drained or time.monotonic() >= deadline:
                            break
                        time.sleep(0.1)
                    if not group_drained:
                        group_signal(signal.SIGKILL)
                    # A zombie-only group cannot fork. Retain the anchor until
                    # the final census; do not signal it merely for being present.
                    deadline = time.monotonic() + 5
                    while True:
                        rows = process_census()
                        members = observe(rows)
                        group_drained = drained(rows, members)
                        if group_drained or time.monotonic() >= deadline:
                            break
                        time.sleep(0.1)
                    if not group_drained:
                        cleanup_errors.append("group did not drain while unreaped leader anchored its PGID")
                except Exception as error:
                    cleanup_errors.append(str(error))
                    try:
                        group_signal(signal.SIGKILL)
                    except Exception as nested:
                        cleanup_errors.append(str(nested))
                try:
                    code = child.wait(timeout=5)
                    reaped = True
                except subprocess.TimeoutExpired:
                    cleanup_errors.append("direct child remained unreaped after group kill")
                try:
                    deadline = time.monotonic() + 2
                    while True:
                        rows = process_census()
                        remaining = {pid: row for (pid, birth), row in seen.items()
                                     if pid in rows and rows[pid]["lstart"] == birth}
                        remaining.update({pid: row for pid, row in rows.items() if row["pgid"] == child.pid})
                        if not remaining or time.monotonic() >= deadline:
                            break
                        time.sleep(0.1)
                    owner["remaining_observed"] = list(remaining.values())
                    retirement_census_complete = True
                except Exception as error:
                    cleanup_errors.append(str(error))
    try:
        after = {str(path): identity(path) for path in inputs}
        write(folder / "inputs-after.json", after)
        unchanged = before == after and observe_source() == source
    except Exception as error:
        unchanged = False
        reason = reason or f"final input/source seal incomplete: {error}"
        write(folder / "inputs-after.json", {"error": str(error)})
    if STOP is not None:
        reason = reason or f"signal {STOP}"
    owner.update(status="terminal", child_exit=code, stop_reason=reason,
                 wall_seconds=time.monotonic() - started, peak_sampled_rss_kib=peak,
                 retirement_census_complete=retirement_census_complete,
                 observed_births=[{**row, "absent": (not any(item["pid"] == row["pid"] and item["lstart"] == row["lstart"]
                                                             for item in owner["remaining_observed"]))
                                                     if retirement_census_complete else None} for row in seen.values()],
                 cleanup_errors=cleanup_errors, whole_inputs_unchanged=unchanged,
                 max_observed_sample_interval_seconds=max([b-a for a,b in zip(times, times[1:])] or [0]),
                 containment_complete=bool(child is not None and reaped and group_drained and retirement_census_complete and not escaped
                                           and not cleanup_errors and not owner["remaining_observed"]),
                 observation_complete=reason is None and unchanged)
    try:
        write(folder / "owner.json", owner)
    finally:
        for number, handler in handlers.items():
            signal.signal(number, handler)
    return owner


def prepare_metal_source(archive, output, modules, remove_forced_inline=False):
    """Recover the exact pair and, when requested, its annotation-only control."""
    require(identity(archive) == {"bytes": 9362078,
            "sha256": "12060f52a9cadba5dabedb5a99b5f4293893bd066d544f2ba4e58c612d2357db"},
            "Metal source archive differs from pinned artifact 11876607522")
    prefix = ("out/macos-retained-diagnostic/observation/"
              "1-PhysicalInitializationRetainsTheHamiltonianResidual/shaders/")
    folder = output / "metal-source"
    folder.mkdir()
    with zipfile.ZipFile(archive) as bundle:
        require(len(bundle.namelist()) == len(set(bundle.namelist())), "duplicate Metal archive member")
        for suffix in ["spv", "metal"]:
            (folder / ("initialize." + suffix)).write_bytes(
                bundle.read(prefix + "shader-cs-94cd92f131b7e6d1." + suffix))
    require(identity(folder / "initialize.spv")["sha256"] ==
            modules["retained_program::kInitializePortableShader"]["sha256"],
            "captured Initialize SPIR-V differs from the original embedded module")
    require(identity(folder / "initialize.metal") == {"bytes": MSL_BYTES, "sha256": MSL_SHA},
            "captured Initialize Metal source differs")
    derived = folder / "initialize-without-forced-inline.metal"
    if remove_forced_inline:
        derived.write_bytes(unforced_msl((folder / "initialize.metal").read_bytes()))
    write(folder / "provenance.json", {"artifact": 11876607522, "archive": identity(archive),
          "source_run": 38526893803, "source_job": 115729265254,
          "prior_bundle": "attestations/native-vulkan/4bcd987/macos-retained-diagnostic-6e13108",
          "prior_evidence_commit": "95ddf4aa1dba0dcec0ef8ca96ea8694d12c82935",
          "spirv": identity(folder / "initialize.spv"), "metal": identity(folder / "initialize.metal"),
          "annotation_only_source_control": {"path": str(derived), **identity(derived),
              "removed_forced_inline_declarations": FORCED_INLINE_COUNT,
              "comparison_control": "Exact captured MSL with the same selected compiler options"} if remove_forced_inline else None,
          "scope": "Exact prior dumped pair; synchronous native Metal diagnostic only, without dispatch or numerical qualification."})
    return folder / "initialize.metal"


def case_result(path, case):
    try:
        tree = ET.parse(path).getroot()
        cases = tree.findall(".//testcase")
        require(int(tree.get("tests", "-1")) == 1 and len(cases) == 1, "expected exactly one GoogleTest case")
        node = cases[0]
        require(node.get("classname") + "." + node.get("name") == case, "different GoogleTest case")
        require(node.get("status") == "run" and node.get("result") == "completed", "case did not complete")
        require(not any(tree.findall(".//" + tag) for tag in ["failure", "error", "skipped"]), "case failed or skipped")
        for parent in [tree, *tree.findall(".//testsuite")]:
            require(all(int(parent.get(key, "0")) == 0 for key in ["failures", "errors", "skipped", "disabled"]), "nonzero GoogleTest counters")
        return {"passed": True}
    except (OSError, ValueError, TypeError, ET.ParseError) as error:
        return {"passed": False, "reason": str(error)}


def dump_join(folder, modules):
    hashes = {}
    for name, pin in modules.items():
        hashes.setdefault(pin["sha256"], []).append(name)
    files = {}
    for path in sorted(folder.rglob("*")):
        if path.is_file():
            pin = identity(path)
            item = dict(pin)
            if path.suffix.lower() == ".spv":
                item["matching_embedded_arrays"] = hashes.get(pin["sha256"], [])
            files[str(path.relative_to(folder))] = item
    spirv = [pin for pin in files.values() if "matching_embedded_arrays" in pin]
    return {"files": files, "all_dumped_spirv_bound": bool(spirv) and all(pin["matching_embedded_arrays"] for pin in spirv),
            "scope": "Byte-exact dumped SPIR-V aliases; raw pipeline-list/MSL associations retained. Loaded modules alone do not prove execution or successful compilation."}


def owner_self_test(output, environment, inputs, *, source_snapshot=None):
    commands = [("exit-zero", "pass", 5, 0, None),
                ("exit-seven", "raise SystemExit(7)", 5, 7, None),
                ("timeout-tree", "import subprocess,sys,time; subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)']); time.sleep(30)", 1, None, "deadline")]
    for label, code, seconds, expected_exit, reason in commands:
        result = run_owned([sys.executable, "-c", code], environment, output / label, inputs, seconds=seconds, source_snapshot=source_snapshot)
        require(result["containment_complete"] and result["whole_inputs_unchanged"], "owner cleanup self-test failed")
        require(result["stop_reason"] == reason and (expected_exit is None or result["child_exit"] == expected_exit), "owner outcome self-test failed")
        if label == "timeout-tree":
            require(len(result["observed_births"]) >= 2, "timeout self-test never observed descendant")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--prepare-only", action="store_true", help="verify/stage the immutable export without native execution")
    parser.add_argument("--metal-source-archive", type=Path,
                        help="diagnose exact captured Initialize MSL through native Metal instead of rerunning controls")
    parser.add_argument("--metal-optimization-size", action="store_true",
                        help="change only native Metal library optimization from Default to Size")
    parser.add_argument("--metal-remove-forced-inline", action="store_true",
                        help="remove only captured forced-inline annotations, retaining the selected compiler options")
    args = parser.parse_args()
    require(not args.metal_optimization_size or args.metal_source_archive,
            "Metal Size diagnostic requires the pinned Metal source archive")
    require(not args.metal_remove_forced_inline or args.metal_source_archive,
            "annotation-only source control requires the pinned archive")
    memory_self_test()
    output = args.output.resolve()
    require(output.is_relative_to(ROOT / "out") and not output.exists(), "output must be a fresh directory under out/")
    source_state()
    require(args.prepare_only or sys.platform == "darwin", "runtime diagnostic must execute natively on macOS")
    output.mkdir(parents=True)
    result = {"status": "incomplete", "qualification_claimed": False, "candidate_source": SOURCE,
              "runner_source": source_state(), "host": {"platform": sys.platform, "machine": platform.machine()},
              "ci": {key: os.environ.get(key) for key in ["GITHUB_RUN_ID", "GITHUB_RUN_ATTEMPT", "GITHUB_JOB", "RUNNER_NAME", "RUNNER_OS", "RUNNER_ARCH"]}, "cases": []}
    write(output / "diagnostic.json", result)
    try:
        stage, app, backend, modules = prepare(args.archive.resolve(), output)
        metal_source = prepare_metal_source(args.metal_source_archive.resolve(), output, modules,
                                           args.metal_remove_forced_inline) if args.metal_source_archive else None
        if args.prepare_only:
            result["status"] = "prepared_without_execution"
            return 0
        environment = dict(os.environ)
        # No inherited case selection, precision override or test randomization.
        environment = {key: value for key, value in environment.items()
                       if not key.startswith(("GTEST_", "SIRIUS_", "MVK_CONFIG_"))
                       and key not in {"VK_ICD_FILENAMES", "VK_ADD_DRIVER_FILES",
                                       "VK_LOADER_DRIVERS_SELECT", "VK_LOADER_DRIVERS_DISABLE"}}
        environment.update(LC_ALL="C", DYLD_PRINT_LIBRARIES="1", MVK_CONFIG_DEBUG="0",
                           MVK_CONFIG_TRACE_VULKAN_CALLS="3")
        icd = Path(environment["VK_DRIVER_FILES"])
        require(icd.is_file(), "explicit MoltenVK ICD is missing")
        require(not any(environment.get(key) for key in ["VK_ICD_FILENAMES", "VK_ADD_DRIVER_FILES"]), "ambiguous provider selector")
        inputs = [args.archive.resolve(), Path(__file__).resolve(), ROOT / ".github/workflows/ci.yml",
                  ROOT / "scripts/runtime_identity.py", ROOT / "scripts/governance/build_policy.py",
                  ROOT / "scripts/macos_metal_probe.py",
                  output / "preparation.json", output / "embedded-modules.json"]
        if metal_source:
            inputs += [args.metal_source_archive.resolve(), *metal_source.parent.iterdir()]
        inputs += [path for path in stage.rglob("*") if path.is_file()]
        library = json.loads(icd.read_text())["ICD"]["library_path"]
        driver = (icd.parent / library).resolve()
        sdk = Path(environment["VULKAN_SDK"])
        loaders = [path.resolve() for path in [sdk / "lib/libvulkan.1.dylib", sdk / "lib/VulkanLoader/lib/libvulkan.1.dylib"] if path.is_file()]
        inputs += [icd.resolve(), driver, *loaders]
        (output / "owner-self-test").mkdir()
        owner_self_test(output / "owner-self-test", environment, inputs)
        def inventory(label):
            receipt = run_owned([str(app), "--json", "info", "system"], environment, output / label, inputs, seconds=60)
            require(receipt["containment_complete"] and receipt["observation_complete"] and receipt["child_exit"] == 0, "native inventory query incomplete")
            data = json.loads((output / label / "stdout.log").read_text())
            require(str(data.get("platform")).casefold() == "macos" and data.get("wsl2") is False, "inventory is not native macOS")
            _, selected = selected_inventory(data, True)
            require(selected["kind"] in ["integrated", "discrete"] and "moltenvk" in (selected["driver_name"] + " " + selected["driver_info"]).lower(), "selected device is not physical MoltenVK")
            return data, selected
        initial, selected = inventory("inventory-initial")
        environment["SIRIUS_VULKAN_DEVICE"] = str(selected["index"])
        bound, selected = inventory("inventory-selected")
        require(initial == bound, "explicit selection changed the observed inventory")
        providers = provider_inputs(app, environment, selected, True)
        write(output / "runtime-inputs.json", {"inventory": bound, "provider_inputs": providers,
                                               "source": result["runner_source"], "candidate_revision": SOURCE,
                                               "qualification_claimed": False})
        inputs += [output / "runtime-inputs.json"]
        if metal_source:
            result["diagnostic_mode"] = "native_metal_pipeline_error"
            folder = output / "metal-pipeline"
            command = [sys.executable, "-B", str(ROOT / "scripts/macos_metal_probe.py"),
                       "--msl", str(metal_source), "--expected-device-name", selected["name"]]
            if args.metal_optimization_size:
                command += ["--optimization-level", "size"]
            if args.metal_remove_forced_inline:
                command += ["--without-forced-inline", str(metal_source.with_name("initialize-without-forced-inline.metal"))]
            receipt = run_owned(command, environment, folder, inputs)
            require(receipt["containment_complete"] and receipt["observation_complete"],
                    "native Metal observation/containment incomplete")
            events = [json.loads(line) for line in (folder / "stdout.log").read_text().splitlines()]
            if args.metal_remove_forced_inline:
                source_events = [event for event in events if event.get("event") == "input"]
                option_events = [event for event in events if event.get("event") == "compile_options"]
                derived = metal_source.with_name("initialize-without-forced-inline.metal")
                require(len(source_events) == 1 and len(option_events) == 1,
                        "annotation control lacks source/options observations")
                observed, options = source_events[0], option_events[0]
                require(observed["sha256"] == MSL_SHA and observed["bytes"] == MSL_BYTES
                        and observed["compiled_path"] == str(derived)
                        and {"bytes": observed["compiled_bytes"], "sha256": observed["compiled_sha256"]} == identity(derived)
                        and observed["captured_source_modified"] is True
                        and observed["removed_forced_inline_declarations"] == FORCED_INLINE_COUNT,
                        "native compiler did not consume the sealed annotation-only source")
                require(options["optimization_level"] == (1 if args.metal_optimization_size else 0)
                        and options["changed_parameter"] == {
                            "forced_inline_declarations": {"baseline": FORCED_INLINE_COUNT, "diagnostic": 0}},
                        "native annotation control changed its comparison options")
            require(events and events[-1].get("event") == "terminal"
                    and events[-1].get("qualification_claimed") is False
                    and events[-1].get("execution_claimed") is False,
                    "native Metal probe did not produce a terminal diagnostic")
            created = events[-1].get("pipeline_created")
            require(type(created) is bool and receipt["child_exit"] == (0 if created else 1),
                    "native Metal terminal verdict/exit differs")
            result.update(status="metal_pipeline_observation_complete", numerics_executed=False,
                          metal_probe={"pipeline_created": created, "events": events,
                                       "route": "synchronous native Metal, distinct from MoltenVK's asynchronous route"})
            return receipt["child_exit"]
        for number, case in enumerate(CASES, 1):
            folder = output / (str(number) + "-" + case)
            local = dict(environment, GTEST_OUTPUT="xml:" + str(folder / "gtest.xml"),
                         MVK_CONFIG_SHADER_DUMP_DIR=str(folder / "shaders"))
            # MoltenVK creates the final dump directory inside the fresh case.
            full_name = "RetainedComputeTest." + case
            receipt = run_owned([str(backend), "--gtest_filter=" + full_name], local, folder, inputs)
            require(STOP is None, "diagnostic interrupted; sequence stopped")
            observed = case_result(folder / "gtest.xml", full_name)
            dumps = dump_join(folder / "shaders", modules)
            write(folder / "dump-joins.json", dumps)
            result["cases"].append({"case": full_name, "owner_complete": receipt["containment_complete"] and receipt["observation_complete"],
                                    "test": observed, "exit": receipt["child_exit"],
                                    "shader_capture_bound": dumps["all_dumped_spirv_bound"]})
            write(output / "diagnostic.json", result)
            require(receipt["containment_complete"] and receipt["observation_complete"], "case observation/containment incomplete; sequence stopped")
            require(dumps["all_dumped_spirv_bound"], "shader dumps missing or not bound to original embedded arrays")
        result["status"] = "observations_complete"
        result["all_four_controls_passed"] = all(row["test"]["passed"] and row["exit"] == 0 for row in result["cases"])
        return 0 if result["all_four_controls_passed"] else 1
    except Exception as error:
        result["error"] = f"{type(error).__name__}: {error}"
        print(result["error"], file=sys.stderr)
        return 1
    finally:
        if STOP is not None:
            result["status"] = "incomplete"
            result["error"] = f"diagnostic interrupted by signal {STOP}"
        result["shader_capture_complete"] = len(result["cases"]) == 4 and all(row["shader_capture_bound"] for row in result["cases"])
        if len(result["cases"]) == 4:
            result["all_four_controls_passed"] = all(row["owner_complete"] and row["test"]["passed"] and row["exit"] == 0 for row in result["cases"])
        write(output / "diagnostic.json", result)
        print(json.dumps(result), flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
