"""Seal one short stall-localization diagnostic and its actual debugger inputs."""
from pathlib import Path
import json, re, subprocess

R = Path.cwd(); W = Path(__file__).resolve().parent; P = W.parent
original = json.loads((P / 'runtime/actions.json').read_text())['control-8']
provider = json.loads((P / 'providers.json').read_text())['providers']['candidate']
inputs = set(original['frozen_inputs'])
inputs |= {str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
inputs |= {'out/native-cooperative-transport-intergroup/runtime/control-8/owner.json',
    'out/native-cooperative-transport-intergroup/runtime/control-8/inputs-before.json',
    'out/native-cooperative-transport-intergroup/runtime/control-8/inputs-after.json', '/usr/bin/gdb'}
output = subprocess.check_output(['ldd', '/usr/bin/gdb'], text=True)
inputs |= {str(Path(x).resolve()) for x in re.findall(r'(/[^\s()]+)', output) if Path(x).is_file()}
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
action = {'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'kind': 'diagnostic', 'expected_exit': 0,
    'limits': {'seconds': 35, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'argv': ['/usr/bin/python3', '-B', str(W / 'capture.py'), '{folder}'],
    'frozen_inputs': sorted(inputs), 'software_provider': True,
    'runtime_icd': provider['ICD'], 'runtime_provider': provider['library'],
    'scope': 'One partial GDB localization only, intentionally interrupted after12s for all-thread stacks and killed. Original candidate failedC8 remains rejected; no acceptance, timing comparison, extra benefit slot, retry or frame verdict.'}
target = W / 'actions.json'; assert not target.exists()
target.write_text(json.dumps({'snapshot': action}, indent=2) + '\n')
print(json.dumps({'partial_stack_snapshot_only': True, 'intentional_interrupt_seconds': 12,
                  'outer_seconds': 35, 'acceptance_runs': 0}))
