#!/usr/bin/env python3
"""One-off harmless CMake output-capture and owned-stop diagnostic."""
import ctypes
import hashlib
import json
import os
from pathlib import Path
import selectors
import shutil
import signal
import subprocess
import time

root = Path(__file__).resolve().parent
cmake = shutil.which("cmake")
assert cmake
version = subprocess.run([cmake, "--version"], capture_output=True, timeout=5, check=True)
completed = subprocess.run([cmake, "-P", str(root / "capture.cmake")],
                           capture_output=True, timeout=5)
(root / "capture-stdout.raw").write_bytes(completed.stdout)
(root / "capture-stderr.raw").write_bytes(completed.stderr)
assert completed.returncode == 0, completed.stderr
assert b"CAPTURE_UNCHANGED" in completed.stdout
assert completed.stdout.count(b"-- stdout: a;b [quoted]\n") == 1
assert completed.stderr.count(b"stderr: c;d [quoted]\n") == 1
assert b"child.cmake" in completed.stdout

# Only this diagnostic process adopts its own stopped CMake descendants.
assert ctypes.CDLL(None).prctl(36, 1, 0, 0, 0) == 0

def stat(pid):
    path = Path("/proc") / str(pid) / "stat"
    try:
        fields = path.read_text().rsplit(")", 1)[1].split()
        return {"pid": pid, "start_ticks": int(fields[19]),
                "pgid": int(fields[2]), "sid": int(fields[3])}
    except (FileNotFoundError, ProcessLookupError):
        return None

def members(pgid):
    result = []
    for path in Path("/proc").glob("[0-9]*/stat"):
        item = stat(int(path.parent.name))
        if item and item["pgid"] == pgid:
            result.append(item)
    return result

started = time.monotonic()
child = subprocess.Popen([cmake, "-P", str(root / "stop.cmake")],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         start_new_session=True)
leader = stat(child.pid)
assert leader and leader["pgid"] == child.pid and leader["sid"] == child.pid
selector = selectors.DefaultSelector()
selector.register(child.stdout, selectors.EVENT_READ, "stdout")
selector.register(child.stderr, selectors.EVENT_READ, "stderr")
streams = {"stdout": bytearray(), "stderr": bytearray()}
early = False
owned_before = []
reaped = []
try:
    deadline = started + 5
    while time.monotonic() < deadline:
        for key, _ in selector.select(0.1):
            data = os.read(key.fileobj.fileno(), 65536)
            if data:
                streams[key.data].extend(data)
            else:
                selector.unregister(key.fileobj)
        early = (b"early-child.cmake" in streams["stdout"] and
                 b"EARLY_STDOUT" in streams["stdout"] and
                 b"EARLY_STDERR" in streams["stderr"])
        if early:
            break
finally:
    owned_before = members(child.pid)
    if child.poll() is None:
        os.killpg(child.pid, signal.SIGTERM)
    try:
        child.wait(timeout=2)
    except subprocess.TimeoutExpired:
        os.killpg(child.pid, signal.SIGKILL)
        child.wait(timeout=2)
    # Popen has reaped the direct leader before group waitpid may run.
    deadline = time.monotonic() + 3
    while time.monotonic() < deadline:
        while True:
            try:
                pid, status = os.waitpid(-child.pid, os.WNOHANG)
            except ChildProcessError:
                break
            if pid == 0:
                break
            reaped.append({"pid": pid, "wait_status": status})
        remaining = members(child.pid)
        if not remaining:
            break
        os.killpg(child.pid, signal.SIGKILL)
        time.sleep(0.01)
    out, err = child.communicate(timeout=2)
    streams["stdout"].extend(out)
    streams["stderr"].extend(err)
    selector.close()

remaining = members(child.pid)
(root / "stop-stdout.raw").write_bytes(streams["stdout"])
(root / "stop-stderr.raw").write_bytes(streams["stderr"])
assert early, "early mirrored command/stdout/stderr was not seen before stop"
assert not remaining, remaining
assert child.returncode < 0
for item in owned_before:
    current = stat(item["pid"])
    assert current is None or current["start_ticks"] != item["start_ticks"]
report = {
    "scope": "Harmless CMake diagnostics only; no Sirius render/build/GPU or permanent test",
    "cmake": cmake, "cmake_version": version.stdout.decode(),
    "capture_case_exit": completed.returncode,
    "original_and_mirrored_captured_stdout_stderr_results_exactly_equal": True,
    "mirrored_stdout_stderr_observed_once": True,
    "command_and_early_stdout_stderr_received_before_intentional_owned_stop": early,
    "owned_leader": leader, "owned_group_before_stop": owned_before,
    "leader_exit": child.returncode, "reaped_adopted_children": reaped,
    "remaining_owned_group_members": remaining,
    "elapsed_stop_diagnostic_seconds": time.monotonic() - started,
    "scope_limit": "Only the recorded diagnostic group/identities; not a Sirius performance or qualification result",
}
(root / "diagnostic-result.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
