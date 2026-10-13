"""Owned final-SPIR-V inspection; original v1 compilation records stay intact."""
import ctypes
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

sys.dont_write_bytecode = True
import screen as original

WORK = original.WORK
REPO = original.REPO
GUARD_SOURCE = REPO / "attestations/software-vulkan/8a27321/vulkan-descriptor-reuse/files/out/vulkan-descriptor-reuse/owned_guard.py"
guard = original.module("established_owned_guard", GUARD_SOURCE)
assert ctypes.CDLL(None, use_errno=True).prctl(36, 1, 0, 0, 0) == 0
BOOT = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
CONTROLLER = guard.proc(os.getpid())


def execute(argv, **kwargs):
    index = len(list((WORK / "actions").glob("*.json")))
    prefix = WORK / "actions" / f"{index:02}"
    started = time.monotonic()
    report = {"argv": list(map(str, argv)), "cwd": str(REPO), "boot_id": BOOT,
              "controller": CONTROLLER, "timeout_seconds": 300,
              "maximum_sampled_rss_kib": 4 * 1024**2,
              "minimum_available_kib": 2 * 1024**2,
              "established_group_helper": original.record(GUARD_SOURCE),
              "runner": original.record(Path(__file__)), "observed_births": [],
              "peak_sampled_rss_kib": 0, "termination_reason": None}
    report["inputs_before"] = {str(Path(value).resolve()): original.record(Path(value))
                               for value in argv if Path(value).is_file()}
    handlers = {number: signal.getsignal(number) for number in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)}
    stop = []
    def interrupted(number, frame):
        stop.append(number)
    child = None
    birth = None
    seen = {}
    reaped = []
    failure = None
    cleanup_errors = []
    try:
        for number in handlers:
            signal.signal(number, interrupted)
        assert guard.available_kib() >= report["minimum_available_kib"]
        with prefix.with_suffix(".stdout").open("wb") as stdout, prefix.with_suffix(".stderr").open("wb") as stderr:
            child = subprocess.Popen(argv, cwd=REPO, stdout=stdout, stderr=stderr,
                                     start_new_session=True)
            birth = guard.proc(child.pid)
            report.update(pid=child.pid, saved_start_ticks=birth["start_ticks"])
            while True:
                members = guard.group_members(child.pid)
                for member in members:
                    seen[(member["pid"], member["start_ticks"])] = member
                rss = sum(member["rss_kib"] for member in members)
                report["peak_sampled_rss_kib"] = max(report["peak_sampled_rss_kib"], rss)
                assert not stop, "supervisor_signal"
                assert time.monotonic() - started <= 300, "elapsed_guard"
                assert rss <= report["maximum_sampled_rss_kib"], "rss_guard"
                assert guard.available_kib() >= report["minimum_available_kib"], "reserve_guard"
                code = child.poll()
                if code is not None:
                    break
                time.sleep(0.01)
            report["return_code"] = code
            assert code == 0, "nonzero_child_exit"
    except BaseException as error:
        failure = error
        report["termination_reason"] = str(error)
    finally:
        try:
            if child is not None:
                try:
                    guard.terminate_group(child, birth["start_ticks"] if birth else None, reaped)
                except BaseException as error:
                    cleanup_errors.append("group cleanup: " + repr(error))
                    if failure is None:
                        failure = error
                    # Even unavailable birth metadata must not abandon the
                    # known unreaped Popen leader. This does not prove or
                    # substitute for descendant closure after a group failure.
                    try:
                        if child.poll() is None:
                            child.kill()
                    except BaseException as kill_error:
                        cleanup_errors.append("known leader kill: " + repr(kill_error))
                finally:
                    # Popen owns the known leader's wait status even if group
                    # termination fails; never infer its status from a PID scan.
                    try:
                        report["return_code"] = child.wait(timeout=10)
                    except BaseException as error:
                        cleanup_errors.append("leader wait: " + repr(error))
                        if failure is None:
                            failure = error
                try:
                    report["remaining_owned_group"] = guard.group_members(child.pid)
                    if report["remaining_owned_group"]:
                        cleanup_errors.append("owned group remains")
                        if failure is None:
                            failure = RuntimeError("owned group remains")
                except BaseException as error:
                    report["remaining_owned_group"] = None
                    cleanup_errors.append("final owned group scan: " + repr(error))
                    if failure is None:
                        failure = error
            try:
                for path, expected in report["inputs_before"].items():
                    assert original.record(Path(path)) == expected
            except BaseException as error:
                if failure is None:
                    failure = error
            if stop and failure is None:
                failure = RuntimeError("supervisor_signal")
            report["termination_reason"] = repr(failure) if failure is not None else None
            report["cleanup_errors"] = cleanup_errors
            report["supervisor_signals"] = stop
            report["observed_births"] = list(seen.values())
            report["reaped_descendants"] = reaped
            report["elapsed_seconds"] = time.monotonic() - started
            report["stdout"] = original.record(prefix.with_suffix(".stdout")) if prefix.with_suffix(".stdout").exists() else None
            report["stderr"] = original.record(prefix.with_suffix(".stderr")) if prefix.with_suffix(".stderr").exists() else None
            prefix.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
            # Keep the stop handlers active through cleanup, postseals and
            # publication, and reject a signal observed during receipt writing.
            if stop and failure is None:
                failure = RuntimeError("supervisor_signal_during_publication")
                report["termination_reason"] = repr(failure)
                prefix.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
        finally:
            for number, handler in handlers.items():
                signal.signal(number, handler)
    if failure is not None:
        raise failure
    print(f"owned action {index}: exit0, group absent, sampled peak {report['peak_sampled_rss_kib']}KiB", flush=True)
    return subprocess.CompletedProcess(argv, report["return_code"])


def inspect():
    for name, expected in json.loads((WORK / "source-freeze.json").read_text()).items():
        assert original.record(WORK / name) == expected, name
    builder = original.module("current_builder", REPO / "scripts/build-retained-kernels.py")
    cache = (REPO / "bin/linux-gcc/CMakeCache.txt").read_text()
    disassembler = re.search(r"^SIRIUS_SPIRV_DIS:FILEPATH=(.*)$", cache, re.M)[1]
    results = json.loads((WORK / "compiler-result.json").read_text())
    destination = WORK / "reinspection"
    destination.mkdir()
    checked = []
    for result in results["modules"]:
        variant = result["variant"]
        normal = variant.endswith("-normal")
        text_by_owner = {}
        for owner in ("stock", "candidate"):
            spv = WORK / "modules" / (owner + "-" + variant + ".spv")
            assert original.record(spv) == result[owner]
            assembly = destination / (owner + "-" + variant + ".spvasm")
            execute([disassembler, str(spv), "-o", str(assembly)])
            text = assembly.read_text()
            text_by_owner[owner] = text
            caps = re.findall(r"OpCapability (\S+)", text)
            assert caps == (["Shader", "RoundingModeRTE"] if normal else ["Shader"])
            assert set(re.findall(r"OpTypeInt (\d+) [01]", text)) == {"32"}
            assert not re.search(r"\b(?:FPFastMathMode|FPFastMathDefault|RelaxedPrecision|OpFmaKHR)\b", text)
            modes = re.findall(r"OpExecutionMode %\S+ (\S+)([^\n]*)", text)
            assert [kind for kind, _ in modes] == (["LocalSize", "RoundingModeRTE"] if normal else ["LocalSize"])
            if normal:
                assert modes[-1] == ("RoundingModeRTE", " 32")
            floating = re.findall(r"(%\S+) = (OpF(?!unction)\S+) (%\S+) ", text)
            if owner == "candidate":
                assert re.findall(r"OpTypeFloat (\d+)", text) == ["32"]
                float_type = re.search(r"(%\S+) = OpTypeFloat 32", text)[1]
                mul = [(value, kind) for value, op, kind in floating if op == "OpFMul"]
                assert len(mul) == 1 and mul[0][1] == float_type
                assert "OpDecorate " + mul[0][0] + " NoContraction" in text
                bodies = re.findall(r"%\S+ = OpFunction .*?OpFunctionEnd", text, re.S)
                body = [body for body in bodies if "= OpFMul " in body]
                assert len(body) == 1 and re.match(r"%PB32ProductResidual(?:_\d+)? = OpFunction ", body[0])
                filtered = re.sub(r"(?m)^.* = OpFMul .*\n", "", text)
                assert filtered.count("OpFMul") == 0
                if not normal:
                    assert len(floating) == 1
            else:
                filtered = text
                if not normal:
                    assert not floating and "OpTypeFloat" not in text
            if normal:
                # Remove only the two already-qualified mode declarations for
                # the production structural checker; this text is never assembled.
                filtered = re.sub(r"(?m)^.*OpCapability RoundingModeRTE\n", "", filtered)
                filtered = re.sub(r"(?m)^.*OpExecutionMode %\S+ RoundingModeRTE 32\n", "", filtered)
                builder.portable_normal_sum_controls(filtered)
        interfaces = lambda text: sorted(re.findall(r"OpDecorate %\S+ (DescriptorSet|Binding|ArrayStride) (\d+)", text))
        assert interfaces(text_by_owner["stock"]) == interfaces(text_by_owner["candidate"])
        checked.append({"variant": variant, "stock_spv": result["stock"], "candidate_spv": result["candidate"],
                        "stock_disassembly": original.record(destination / ("stock-" + variant + ".spvasm")),
                        "candidate_disassembly": original.record(destination / ("candidate-" + variant + ".spvasm")),
                        "final_types_modes_product_fmul_and_original_float_bodies_checked": True,
                        "binding_stride_decoration_inventory_exact": True})
    for name, expected in json.loads((WORK / "source-freeze.json").read_text()).items():
        assert original.record(WORK / name) == expected, name
    original.write("final-spv-inspection.json", {"passed": True,
        "scope": "Fresh final-binary disassembly and bounded structural checks. Guard/correction semantics require independent compiled review; binding inventories alone are not complete interface equivalence.",
        "original_v1_owners": {"count": 60, "normal_exit_only": True,
             "sampled_rss_reserve_and_descendant_closure_not_recorded": True},
        "compiled_original_runner": original.record(WORK / "screen.py"),
        "inspection_runner": original.record(Path(__file__)), "modules": checked})


if __name__ == "__main__":
    inspect()
