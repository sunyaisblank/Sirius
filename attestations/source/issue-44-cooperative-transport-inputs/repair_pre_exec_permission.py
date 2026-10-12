"""Preserve a zero-child launch failure and restore archived executable mode."""
from pathlib import Path
import hashlib
import json
import os
import stat

root = Path.cwd()
work = Path(__file__).resolve().parent / 'benefit'
identity = lambda p: {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
owner_path = work / 'a1/owner.json'
owner = json.loads(owner_path.read_text())
sequence = json.loads((work / 'software-sequence.json').read_text())
freeze = json.loads((work / 'retention-gate-frozen.json').read_text())
assert owner['status'] == 'terminal' and owner['child_exit'] is None
assert owner['stop_reason'] == 'observer_exception'
assert owner['observer_error'] == "PermissionError(13, 'Permission denied')"
assert not owner.get('child_pid') and not owner['observed_births']
assert not owner['provider_mappings'] and not owner['remaining_owned_processes']
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['cleanup_errors']
assert (work / 'a1/samples.jsonl').stat().st_size == 0
assert not (work / 'a1/gtest.xml').exists()
assert sequence['results'] == [] and not sequence['completed']
assert not any((work / label).exists() for label in ('b1', 'b2', 'a2'))
for producer, value in freeze['executables'].items():
    assert identity(root / value['path']) == {key: value[key] for key in ('bytes', 'sha256')}
baseline = root / freeze['executables']['baseline']['path']
assert stat.S_IMODE(baseline.stat().st_mode) == 0o644
assert not os.access(baseline, os.X_OK)
frozen = [work / p for p in ('retention-gate-frozen.json', 'actions.json', 'run_sequence.py', 'compare_software.py', 'software_contracts.py')]
before = {str(p.relative_to(root)): identity(p) for p in frozen}
failed = [*sorted((work / 'a1').rglob('*')), work / 'software-sequence.json']
manifest = {str(p.relative_to(work)): identity(p) for p in failed if p.is_file()}
saved = work / 'pre-exec-permission-failure'
assert not saved.exists()
saved.mkdir()
(work / 'a1').rename(saved / 'a1')
(work / 'software-sequence.json').rename(saved / 'software-sequence.json')
assert all(identity(saved / name) == value for name, value in manifest.items())
baseline.chmod(0o755)
assert os.access(baseline, os.X_OK)
value = freeze['executables']['baseline']
assert identity(baseline) == {key: value[key] for key in ('bytes', 'sha256')}
assert before == {str(p.relative_to(root)): identity(p) for p in frozen}
record = {'passed': True, 'preserved_failed_launch_files': manifest,
          'failure_owner_sha256': manifest['a1/owner.json']['sha256'],
          'actual_child_ever_started': False, 'completed_timing_slots': 0,
          'baseline_permission_before': '0644', 'baseline_permission_after': '0755',
          'baseline_bytes_and_sha256_unchanged': True,
          'frozen_inputs_unchanged': before,
          'independent_review': 'command_owner_review confirmed zero-child Popen failure and launch repair, not a performance retry; no adoption clearance.',
          'scope': 'The first actual frozen ABBA may now launch. No timing or performance result was discarded, no source, executable byte, threshold, budget, input or sequence change.'}
(work / 'pre-exec-permission-correction.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'passed': True, 'preserved_failed_files': len(manifest), 'actual_timing_slots_before_repair': 0, 'frozen_inputs_unchanged': True}))
