"""Seal I with unchanged G source,32-bit offsets and OS-guarded AVX2 packets."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-simd-rows/prepare_ispc_compile.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '6df49b3e0481789fc77c4ef04606303d26941b20d002511f5fc74996ba62dc8b'
source = basis.read_text()
address_model = {'pointer_bits': 64, 'offset_bits': 32, 'maximum_capacity': 65536,
                 'maximum_admitted_view_bytes': {}}
for stage, input_stride, output_stride, program_words in (
    ('transport', 230, 4699, 15967), ('endpoint', 225, 3829, 21143),
    ('dopri_phase', 1810, 4149, 10689)):
    extents = [(1 + 65536 * input_stride + program_words) * 4,
               65536 * output_stride * 4]
    assert all(0 < value < 2**31 for value in extents)
    address_model['maximum_admitted_view_bytes'][stage] = extents
preserved = '''inputs = dict(scalar['source_inputs'])
# Bind original G source/objects and the accepted J derivation before K.
for receipt_name in ('ispc-uniform-compile-result.json', 'adapter-compile-result.json'):
    receipt_path = root / 'out/native-cpu-simd-rows' / receipt_name
    previous = read(receipt_path)
    assert previous['all_prepared_source_bytes_unchanged']
    inputs.update(previous['source_inputs'])
    inputs[str(receipt_path.relative_to(root))] = pin(receipt_path)
    for item in previous['objects']:
        assert item['exit'] == 0
        value = item['object']
        inputs[value['path']] = {k: value[k] for k in ('bytes', 'sha256')}
'''
preserved += '''
prior = root / 'out/native-cpu-simd-output-columns'
prior_receipt = prior / 'ispc-uniform-compile-result.json'
previous = read(prior_receipt)
assert previous['all_prepared_source_bytes_unchanged'] and previous['all_objects_compiled']
inputs.update(previous['source_inputs'])
inputs[str(prior_receipt.relative_to(root))] = pin(prior_receipt)
for item in previous['objects']:
    value = item['object']
    inputs[value['path']] = {k: value[k] for k in ('bytes', 'sha256')}
'''
changes = [
    ("derived = read(work / 'ispc-uniform-derivation.json')", "derived = read(work / 'ispc-derivation.json')"),
    ("'ispc-uniform-derivation.json', 'prepare_ispc_compile.py'", "'ispc-derivation.json', 'prepare_ispc_compile.py'"),
    ("destination = root / 'bin/linux-gcc/native-cpu-simd-rows/objects-uniform'",
     "destination = root / 'bin/linux-gcc/native-cpu-simd-structural-uniform/objects-uniform'"),
    ("'--addressing=64'", "'--addressing=32'"),
    ("'--target=sse2-i32x8'", "'--target=avx2-i32x8'"),
    ("'--cpu=x86-64'", "'--cpu=core-avx2'"),
    ("inputs = dict(scalar['source_inputs'])\n", preserved),
    ("old = read(work / 'generation/actions.json')['c-source']",
     "old = read(root / 'out/native-cpu-simd-rows/generation/actions.json')['c-source']"),
    ("actions = read(work / 'compile/actions.json')",
     "assert not (work / 'compile/actions.json').exists()\nactions = {}"),
    ("prepared = {'revision': revision", "prepared = {'address_model': address_model, 'revision': revision"),
    ('Three complete per-row ISPC objects with uniform canonical operands/private scratch indices;',
     'Uniform fixed structural indices and contiguous private output/coefficient accesses, preserving J arithmetic/all9380 calls; required core-avx2 device guarded by OS-enabled AVX2 plus x86-64-v3, with original scalar fallback and unchanged strict FP/math flags;'),
]
reuse_plan = '''prior_derivation = read(root / 'out/native-cpu-simd-output-columns/ispc-derivation.json')
prior_compile = read(root / 'out/native-cpu-simd-output-columns/ispc-uniform-compile-result.json')
prior_objects = {row['stage']: row for row in prior_compile['objects']}
for item in objects:
    slug = item['stage']
    prior_source = prior_derivation['stages'][slug]
    if derived['stages'][slug]['derived'] == prior_source['derived']:
        assert (root / derived['stages'][slug]['derived_path']).read_bytes() == (root / prior_source['derived_path']).read_bytes()
        old = prior_objects[slug]
        assert old['exit'] == 0 and old['argv'][2:-3] == flags
        item['reuse'] = old['object']
        assert pin(root / item['reuse']['path']) == {key: item['reuse'][key] for key in ('bytes', 'sha256')}
'''
changes.append(("for name in ('compiler-acquisition.json', 'compiler-environment.json', 'derive_ispc.py',", reuse_plan + "for name in ('compiler-acquisition.json', 'compiler-environment.json', 'derive_ispc.py',"))
for before, after in changes:
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
