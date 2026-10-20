import ast
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PUB = Path(__file__).resolve().parent
SOURCE = ROOT / 'out/native-coefficient-sharing-runtime-source'
phase, = sys.argv[1:]
assert phase in ('timestamp', 'abba', 'remaining')
WORK = ROOT / ('out/native-coefficient-sharing-runtime-' + phase)


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


plan, binding = read(WORK / 'plan.json'), read(WORK / 'execution-bindings.json')
assert plan['admitted'] is binding['admitted'] is False
assert plan['source_revision'] == binding['source_revision'] == 'e1c618f83861c478fe18c5a22f5d9aa27f6280fd'
assert plan['runner_source_revision'] == binding['runner_source_revision'] == '2c8ef08693d46f57d5dea14efdbddd4874dac684'
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT).decode().strip() == plan['runner_source_revision']
assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT)
assert plan['control_profile'] == binding['control_profile'] == phase
assert (plan['owner_rss_bytes'], plan['minimum_wsl_available_kib'], plan['terminal_observer_seconds']) == (4294967296, 2097152, 60)
expected_order = {'timestamp': ['control-8'], 'abba': ['a1', 'b1', 'b2', 'a2'],
                  'remaining': ['control-' + str(i) for i in range(1, 8)]}[phase]
assert plan['order'] == expected_order and set(binding['actions']) == set(binding['terminal_actions']) == set(expected_order)
for action in expected_order:
    value = plan['actions'][action]
    assert (value['owner_seconds'], value['bridge_outer_seconds'], value['expected_exit']) == (300, 1380, 0)
    expected_arm = 'baseline' if action in ('a1', 'a2') else 'candidate'
    assert value['arm'] == expected_arm
    command = binding['actions'][action]['argv']
    assert command[:3] == ['/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe', '-NoProfile', '-NonInteractive']
    assert value['test'] in command[-1] and 'original_' + expected_arm + '_controls.exe' in command[-1]
assert binding['artifact']['id'] == 12250389795 and binding['artifact']['run_id'] == '39393948651'
assert binding['native_build_gate_authority_used'] is False
for name, row in binding['files'].items():
    path = Path(name) if name.startswith('/') else ROOT / name
    assert pin(path) == {key: row[key] for key in ('bytes', 'sha256')}, name
for name in ['run_system.py', 'run_terminal.py', 'run_native_control.ps1', 'terminal-audit.ps1', 'emergency_native_cleanup.ps1']:
    assert (WORK / name).read_bytes() == (SOURCE / name).read_bytes()
for name in ['native_evidence.py', 'numeric_gates.py', 'compare_native.py', 'read_controls_native.py',
             'expected-readback-layout.json', 'emergency_native_cleanup.ps1', 'run_native_control.ps1', 'terminal-audit.ps1']:
    assert (SOURCE / name).read_bytes() == (SOURCE / 'originals' / name).read_bytes()
for name in ['root-actual-producer-eligibility.json', 'root-publication-check.json', 'root-manual-admission.json']:
    value = read(PUB / name)
    assert value['passed'] and value['runtime_admitted'] is False
if phase != 'timestamp':
    assert read(ROOT / 'out/native-coefficient-sharing-timestamp-evidence/timestamp-acceptance.json')['passed']
if phase == 'remaining':
    assert read(ROOT / 'out/native-coefficient-sharing-ABBA-evidence/native-benefit.json')['passed']
before = {}
after = {}
for name, value in [('plan.json', plan), ('execution-bindings.json', binding)]:
    path = WORK / name
    before[name] = pin(path)
    held = PUB / (phase + '-held-' + name)
    assert not held.exists()
    shutil.copy2(path, held)
    value['admitted'] = True
    path.write_text(json.dumps(value, indent=2) + '\n')
    after[name] = pin(path)
report = {'schema': 'sirius-native-coefficient-sharing-concrete-phase-root-release-v1',
          'passed': True, 'phase': phase, 'original_order': expected_order,
          'input_files_rehashed_equal': len(binding['files']), 'before': before, 'after': after,
          'source_gate': pin(Path(__file__).resolve()),
          'scope': 'Only this concrete original phase/terminals released; no retry, frame, remaining-stage, adoption or qualification right implied.'}
output = PUB / (phase + '-root-release.json')
assert not output.exists()
output.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report))
