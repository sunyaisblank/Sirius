"""Prepare three strict object recompiles using unchanged generated sources."""
from pathlib import Path
import hashlib
import json
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
basis = root / 'out/native-cpu-slang-bridge'
destination = root / 'bin/linux-gcc/native-cpu-active-prefix/typed-batch'


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not (work / 'prepared-batch-inputs.json').exists()
original = json.loads((basis / 'batch-compile-result.json').read_text())
for name, expected in original['source_inputs'].items():
    assert pin(root / name) == expected, name
destination.mkdir(parents=True, exist_ok=False)
inputs = dict(original['source_inputs'])
inputs[str((work / 'native_packet.hpp').relative_to(root))] = pin(work / 'native_packet.hpp')
objects = []
for item in original['objects']:
    assert item['exit'] == 0
    argv = list(item['argv'])
    include = argv.index('-I')
    assert argv[include + 1] == str(basis)
    argv[include:include + 2] = ['-I', str(work), '-I', str(basis)]
    output = destination / Path(item['object']['path']).name
    argv[argv.index('-o') + 1] = str(output)
    objects.append({'stage': item['stage'], 'argv': argv,
                    'output': str(output.relative_to(root))})
prepared = {'inputs': inputs, 'objects': objects,
    'scope': 'Exact original CPP/nativeHPP/canonicalprogramHPP reused. Only native_packet.hpp differs: immutablecapacity defines extents/tail; activeprefix defines dispatchedgroups. No Slang/SPV/host regeneration and no compiler execution by this preparer.'}
(work / 'prepared-batch-inputs.json').write_text(json.dumps(prepared, indent=2) + '\n')
old = json.loads((basis / 'build/actions.json').read_text())['typed-three-stage-batch']
frozen = {name for name in old['frozen_inputs'] if name.startswith('/')}
frozen |= set(inputs)
frozen |= {str(p.relative_to(root)) for p in work.iterdir() if p.is_file()}
frozen |= {str(basis / 'build_batch.py'),
           str(root / 'out/native-cpu-row-ownership/cold-frame/run_owned.py'),
           str(root / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py')}
assert all((root / name).is_file() for name in frozen)
action = {**old,
    'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'argv': ['/usr/bin/python3', '-B', str(work / 'build_batch.py')],
    'frozen_inputs': sorted(frozen),
    'scope': 'One original strictC++26 three-object batch reusing all generated shader/table sources; only activeprefix admission changes. Same baselineISA/SSE and strictFP. No native/provider execution, numerical/benefit/frame acceptance.'}
(work / 'build/actions.json').write_text(json.dumps({'compile': action}, indent=2) + '\n')
print(json.dumps({'prepared_objects': 3, 'sealed_inputs': len(frozen), 'compiled': False}))
