"""Link the new S objects to existing verified diagnostic objects; read exact code."""
from pathlib import Path
import json
import re
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
OLD = R / 'out/native-cpu-accumulate-compatibility'
sys.path.insert(0, str(R / 'out/native-cpu-field-compatibility'))
from elf_proof import pin, associate, rename_proof

mode, folder = sys.argv[1:]
assert mode == 'link' and Path(folder) == W / 'link'
prepared = json.loads((W / 'prepared.json').read_text())
for name, expected in prepared['source_inputs'].items():
    assert pin(R / name) == expected, name
executable = R / prepared['executable']
assert not executable.exists() and not (W / 'link-result.json').exists()
subprocess.run(prepared['link_argv'], cwd=R, check=True)
baseline = R / prepared['baseline']; old = json.loads((OLD / 'compile-completion.json').read_text())
assert rename_proof(R / prepared['original_baseline_object']['path'], baseline,
                    json.loads((OLD / 'prepared.json').read_text())['rename_mapping']) == json.loads((OLD / 'rename-proof.json').read_text())
proof = associate(executable, [{'path': str(baseline), **pin(baseline)},
                              {**prepared['candidate_transport_object'], 'path': str(R / prepared['candidate_transport_object']['path'])}])
def identity(path):
    return {'path': str(path.relative_to(R)), **pin(path)}
proof['executable'] = identity(executable)
proof['objects'][0]['object'] = old['renamed_baseline_object']
proof['objects'][1]['object'] = prepared['candidate_transport_object']
(W / 'link-proof.json').write_text(json.dumps(proof, indent=2) + '\n')
resolution = subprocess.check_output(['/usr/bin/ldd', '-r', str(executable)], text=True, stderr=subprocess.STDOUT)
(W / 'link/runtime-resolution.txt').write_text(resolution)
assert 'undefined symbol' not in resolution and 'not found' not in resolution
libraries = {str(Path(name).resolve(strict=True)) for name in re.findall(r'(/[^\s()]+)', resolution) if Path(name).is_file()}
for name, expected in prepared['source_inputs'].items():
    assert pin(R / name) == expected, name
result = {'passed': True, 'prepared': pin(W / 'prepared.json'),
          'commands': [{'argv': prepared['link_argv'], 'exit': 0}], 'new_CPP_compilations': 0, 'new_links': 1,
          'all8_native_objects': prepared['objects'], 'candidate_transport_object': prepared['candidate_transport_object'],
          'driver_object': old['driver_object'], 'renamed_baseline_object': old['renamed_baseline_object'],
          'rename_proof': old['rename_proof'], 'link_proof': identity(W / 'link-proof.json'),
          'executable': identity(executable), 'runtime_libraries': {name: pin(Path(name)) for name in sorted(libraries)},
          'reused_R_build': prepared['reused_R_build'], 'source_inputs_unchanged': True,
          'scope': prepared['scope'] + ' Every retained function in old/new Transport is associated by relocation-masked opcodes in this exact standalone ELF. Other7 objects pinned/linked; current provider independently joins all8.'}
with (W / 'link-result.json').open('x') as stream:
    json.dump(result, stream, indent=2); stream.write('\n')
print(json.dumps({'linked': True, 'new_CPP_compilations': 0, 'new_links': 1}))
