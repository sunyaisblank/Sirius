import ast
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys
import types

sys.dont_write_bytecode = True
root = Path.cwd()
source = 'src/sirius/kernels/retained_program.py'
baseline_revision = '4d9dfff9390fb65a694fb840821618017c6be4a9'
baseline = types.ModuleType('baseline_program')
candidate = types.ModuleType('candidate_program')
exec(subprocess.check_output(['git', 'show', baseline_revision + ':' + source], text=True), baseline.__dict__)
exec((root / source).read_text(), candidate.__dict__)
decoder_source = root / 'attestations/software-vulkan/b5cab2e/endpoint-schedule-baseline/controls/decode_trial.py'
tree = ast.parse(decoder_source.read_text())
decoder = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'decode')
exec(compile(ast.Module(body=[decoder], type_ignores=[]), str(decoder_source), 'exec'))

def encoded(program, camera=False):
    header = [program['instructions'], program['registers']]
    if not camera:
        header.append(len(program['outputs']))
    return header + program['outputs'] + program['operations'] + program['layer_offsets']

def digest(words):
    return hashlib.sha256(struct.pack('<' + 'I' * len(words), *words)).hexdigest()

unchanged = []
for name in ['camera', 'transport', 'endpoint', 'dense', 'initialize', 'ray_camera', 'dopri_phase', 'schwarzschild_transport']:
    builder = 'build_' + name + '_program'
    a, b = (getattr(m, builder)(True) for m in [baseline, candidate])
    assert a == b, name
    unchanged.append({'program': name, 'encoded_words': len(encoded(a, name in ['camera', 'ray_camera'])),
                      'sha256': digest(encoded(a, name in ['camera', 'ray_camera']))})

reports = []
for module in [baseline, candidate]:
    original = module.compile_program
    captured = []
    def capture(*args, **kwargs):
        result = original(*args, **kwargs)
        captured.append(result['registers'])
        return result
    module.compile_program = capture
    program = module.build_schwarzschild_endpoint_program(True)
    module.compile_program = original
    decoded = decode(program)
    decoded['raw_registers'] = captured[-1]
    decoded['divides'] = program['operations'][::5].count(5)
    decoded['multiplies'] = program['operations'][::5].count(4)
    decoded['encoded_sha256'] = digest(encoded(program))
    decoded['prefix_offsets'] = program['layer_offsets'][:program['layer_offsets'].index(program['prefix_instructions']) + 1]
    reports.append(decoded)
a, b = reports
assert a['root_hashes'][:20] == b['root_hashes'][:20]
assert a['prefix_offsets'] == b['prefix_offsets']
assert a['prefix'] == b['prefix'] == 202
assert (a['instructions'], b['instructions'], a['layers'], b['layers'], a['divides'], b['divides']) == (1943, 1834, 84, 75, 69, 55)
assert a['registers'] == b['registers'] == 612
assert a['outputs'] == b['outputs'] == 100
word_delta = b['words'] - a['words']
assert word_delta == (1834 - 1943) * 5 + (75 - 84) == -554
byte_delta = word_delta * 4
receipt = {'baseline_revision': baseline_revision,
           'candidate_source_sha256': hashlib.sha256((root/source).read_bytes()).hexdigest(),
           'unchanged_encoded_programs': unchanged,
           'schwarzschild_endpoint': {'baseline': a, 'candidate': b},
           'all20_metric_tangent_expression_roots_and_prefix_offsets_exact': True,
           'encoded_prefix_register_identity_not_claimed': True,
           'row_words_unchanged': 3829,
           'immutable_endpoint_table_bytes_delta': byte_delta,
           'independent_fixed_layout_expected_totals': {str(cap): total + byte_delta for cap, total in [(1,704192),(24,3334012),(64,7907612)]},
           'capacity24_endpoint_input_expected_bytes': 106176 + byte_delta,
           'scope': 'Read-only generation and independent serialized DAG/layout checks; no numerical or speedup claim.'}
(Path(__file__).parent/'program-review.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({'unchanged_programs': len(unchanged), 'schwarzschild_instructions': [a['instructions'],b['instructions']],
                  'divides': [a['divides'],b['divides']], 'layers': [a['layers'],b['layers']], 'immutable_delta_bytes': byte_delta,
                  'expected_totals': receipt['independent_fixed_layout_expected_totals'],
                  'endpoint_input': receipt['capacity24_endpoint_input_expected_bytes']}))
