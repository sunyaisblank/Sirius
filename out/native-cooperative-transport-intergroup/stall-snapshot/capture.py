"""One deliberately interrupted debugger snapshot; never an acceptance rerun."""
from pathlib import Path
import hashlib, json, os, signal, subprocess, sys

R = Path.cwd(); W = Path(__file__).resolve().parent; folder = Path(sys.argv[1])
P = W.parent
matched = json.loads((P / 'matched-result.json').read_text())
exe = R / matched['executable']
assert hashlib.file_digest(exe.open('rb'), 'sha256').hexdigest() == matched['executable_sha256']
failed = json.loads((P / 'runtime/control-8/owner.json').read_text())
assert not failed['passed'] and failed['stop_reason'] == 'diagnostic_timeout'
argv = ['/usr/bin/gdb', '-q', '-batch', '-nx',
    '-ex', 'set pagination off', '-ex', 'set confirm off',
    '-ex', 'set print thread-events off', '-ex', 'set debuginfod enabled off',
    '-ex', 'set disable-randomization off', '-ex', 'run',
    '-ex', r'echo \n=== INTENTIONAL PARTIAL DIAGNOSTIC STACK ===\n',
    '-ex', 'thread apply all bt 16', '-ex', 'kill', '--args', str(exe),
    '--gtest_filter=RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults',
    '--gtest_color=no']
interrupted = False
with (folder / 'gdb.log').open('wb') as log:
    child = subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT)
    try:
        try:
            code = child.wait(timeout=12)
        except subprocess.TimeoutExpired:
            interrupted = True
            os.kill(child.pid, signal.SIGINT)
            code = child.wait(timeout=12)
    finally:
        if child.poll() is None:
            child.kill(); child.wait()
text = (folder / 'gdb.log').read_text(errors='replace')
assert interrupted and '#0 ' in text and 'INTENTIONAL PARTIAL DIAGNOSTIC STACK' in text
result = {'snapshot_observed': True, 'intentional_debugger_interrupt_after_seconds': 12,
    'debugger_exit': code, 'executable_sha256': matched['executable_sha256'],
    'failed_original_source': failed['source'],
    'scope': 'One partial debugger localization of a newly observed original C8 timeout. Intentionally interrupted after12s and then killed; no C8, numerical, performance, frame or qualification verdict. Does not replace or retry the failed300s acceptance run.'}
(folder / 'snapshot.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
