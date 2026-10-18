import hashlib
import json
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
config = json.loads((WORK / 'transport-default-config.json').read_text())
unc = config['query_script'].rsplit('\\', 1)[0]
argv = ['/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe',
        '-NoProfile', '-NonInteractive', '-ExecutionPolicy', 'Bypass',
        '-File', unc + '\\guard_query.ps1', '-Config',
        unc + '\\transport-default-config.json']
assert not (WORK / 'transport-default-query').exists()
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip() == config['source_revision']
assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True)
config_hash = hashlib.sha256((WORK / 'transport-default-config.json').read_bytes()).hexdigest()
started = time.monotonic()
with (WORK / 'transport-default-bridge.stdout').open('wb') as out, (WORK / 'transport-default-bridge.stderr').open('wb') as err:
    # The Windows owner independently bounds its child to 120s and 4GiB sampled RSS.
    process = subprocess.Popen(argv, cwd=ROOT, stdout=out, stderr=err)
    birth_ticks = Path('/proc/%d/stat' % process.pid).read_text().split(') ', 1)[1].split()[19]
    result = process.wait()
elapsed = time.monotonic() - started
record = {'argv': argv, 'returncode': result, 'linux_bridge_pid': process.pid,
          'linux_bridge_birth_ticks': birth_ticks, 'linux_bridge_elapsed_seconds': elapsed,
          'config_sha256': config_hash, 'source_revision': config['source_revision'],
          'scope': 'One default Transport live-driver compilation/query; no dispatch, benchmark or numerical acceptance'}
assert hashlib.sha256((WORK / 'transport-default-config.json').read_bytes()).hexdigest() == config_hash
(WORK / 'transport-default-bridge.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record))
raise SystemExit(result)
