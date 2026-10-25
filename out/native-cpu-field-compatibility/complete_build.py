"""Complete code association after the preserved successful compile/link; no build."""
from pathlib import Path
import json
import re
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
from elf_proof import pin, associate, rename_proof

assert sys.argv[1:] == ['proof', str(W / 'proof')]
prepared = json.loads((W / 'prepared.json').read_text())
partial = json.loads((W / 'build-result.json').read_text())
owner = json.loads((W / 'build/owner.json').read_text())
assert not partial['passed'] and len(partial['commands']) == 2
assert all(item['exit'] == 0 for item in partial['commands'])
assert owner['status'] == 'terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['passed'] and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
assert not (W / 'build-completion.json').exists()
for name, expected in prepared['source_inputs'].items():
    assert pin(R / name) == expected, name


def identity(path):
    return {'path': str(path.relative_to(R)), **pin(path)}


OLD = R / 'out/native-cpu-accumulate-compatibility'
old = json.loads((OLD / 'compile-completion.json').read_text())
baseline = R / prepared['baseline']; executable = R / prepared['executable']
assert rename_proof(R / prepared['original_baseline_object']['path'], baseline,
    json.loads((OLD / 'prepared.json').read_text())['rename_mapping']) == json.loads((OLD / 'rename-proof.json').read_text())
associated = [identity(baseline), *prepared['objects'][:2], identity(executable.parent / 'shim.o')]
proof = associate(executable, [{**row, 'path': str(R / row['path'])} for row in associated])
proof['executable'] = identity(executable)
for row, original in zip(proof['objects'], associated):
    row['object'] = original
proof['scope'] = 'All retained functions in original renamed baseline, actual new Transport CPP/layer and diagnostic shim, including unique local-only initialization/cold sections. Declared4B PC relocations masked; relaxation disabled; no numerical claim.'
(W / 'link-proof.json').write_text(json.dumps(proof, indent=2) + '\n')
resolution = subprocess.check_output(['/usr/bin/ldd', '-r', str(executable)], text=True, stderr=subprocess.STDOUT)
(W / 'proof/runtime-resolution.txt').write_text(resolution)
assert 'undefined symbol' not in resolution and 'not found' not in resolution
libraries = {str(Path(name).resolve(strict=True)) for name in re.findall(r'(/[^\s()]+)', resolution) if Path(name).is_file()}
result = {**partial, 'passed': True, 'diagnostic_shim_compilations': 1, 'new_links': 1,
    'additional_compilations_in_proof_completion': 0, 'additional_links_in_proof_completion': 0,
    'executable': identity(executable), 'link_proof': identity(W / 'link-proof.json'),
    'driver_object': old['driver_object'], 'renamed_baseline_object': old['renamed_baseline_object'],
    'rename_proof': old['rename_proof'], 'shim_object': identity(executable.parent / 'shim.o'),
    'candidate_transport_object': prepared['candidate_transport_object'], 'source_inputs_unchanged': True,
    'prior_failed_owner': identity(W / 'build/owner.json'), 'partial_build_receipt': identity(W / 'build-result.json'),
    'runtime_libraries': {name: pin(Path(name)) for name in sorted(libraries)}}
for name, expected in prepared['source_inputs'].items():
    assert pin(R / name) == expected, name
with (W / 'build-completion.json').open('x') as stream:
    json.dump(result, stream, indent=2); stream.write('\n')
print(json.dumps({'proof_complete': True, 'additional_compilations': 0, 'additional_links': 0,
                  'associated_functions': sum(len(row['functions']) for row in proof['objects'])}))
