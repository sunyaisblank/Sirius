"""Read all completed original-builder products; never compile or load a pipeline."""
from pathlib import Path
import hashlib
import json
import re
import struct
import subprocess
import sys
sys.dont_write_bytecode = True
R = Path.cwd().resolve()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/computed-carry-exit-original-prerequisites'
C = R / 'bin/linux-gcc/src/sirius/backend/retained'


def identity(path):
    with Path(path).open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': Path(path).stat().st_size, 'sha256': digest}


def arrays(path):
    result = {}
    for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', path.read_text(), re.S):
        words = [int(v) for v in re.findall(r'(\d+)u', body)]
        assert len(words) == int(count) and name not in result
        result[name] = struct.pack('<' + str(len(words)) + 'I', *words)
    return result


for name, pin in json.loads((W / 'source-freeze.json').read_text()).items():
    assert identity(name) == pin, name
owners = {}
for mode in ['stock', 'candidate']:
    folder = W / ('compile-products-' + mode)
    owner = json.loads((folder / 'owner.json').read_text())
    assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and owner['cleanup_errors'] == [] and owner['remaining_owned_processes'] == []
    assert owner['observed_births'] and all(b['absent'] for b in owner['observed_births'])
    assert json.loads((folder / 'inputs-before.json').read_text()) == json.loads((folder / 'inputs-after.json').read_text())
    owners[mode] = identity(folder / 'owner.json')

expected_changed = {'retained_' + stage + '_portable' + suffix + '.spv'
                    for stage in ['transport', 'endpoint', 'dense', 'initialize', 'dopri_phase']
                    for suffix in ['', '_fp64']}
expected_changed |= {'retained_' + stage + '_portable_normal_sum.spv'
                     for stage in ['transport', 'endpoint', 'dopri_phase']}
expected_changed |= {'retained_' + stage + '_portable_stage7_owner.spv'
                     for stage in ['transport', 'endpoint']}
canonical = {p.name for p in C.glob('*.spv')}
assert len(canonical) == 35 and len(expected_changed) == 15 and expected_changed <= canonical
for mode in ['stock', 'candidate']:
    assert {p.name for p in (B / mode).glob('*.spv')} == canonical
cache = (R / 'bin/linux-gcc/CMakeCache.txt').read_text()
disassembler = re.search(r'^SIRIUS_SPIRV_DIS:FILEPATH=(.*)$', cache, re.M)[1]
shape = lambda text: sorted(re.findall(r'OpExecutionMode %\S+ (.*)|OpDecorate %\S+ ((?:Binding|DescriptorSet|ArrayStride) \d+)|OpMemberDecorate %\S+ (\d+ Offset \d+)', text))
rows = []
for name in sorted(canonical):
    stock, candidate = B / 'stock' / name, B / 'candidate' / name
    assert stock.read_bytes() == (C / name).read_bytes(), name
    changed = stock.read_bytes() != candidate.read_bytes()
    assert changed == (name in expected_changed), name
    texts = [subprocess.check_output([disassembler, str(p)], text=True) for p in [stock, candidate]]
    assert shape(texts[0]) == shape(texts[1]), name
    for opcode in ['OpControlBarrier', 'OpMemoryBarrier']:
        assert texts[0].count(opcode) == texts[1].count(opcode), (name, opcode)
    rows.append({'module': name, 'changed': changed, 'stock': identity(stock), 'candidate': identity(candidate),
                 'stock_equals_canonical': True, 'mode_binding_stride_offset_and_barrier_inventory_equal': True})

headers = {}
for name, count in [('retained_kernels.h', 40), ('retained_owner_kernels.h', 2)]:
    baseline, stock, candidate = [arrays(p / name) for p in [C, B / 'stock', B / 'candidate']]
    assert len(baseline) == count and baseline == stock and baseline.keys() == candidate.keys()
    changed = [key for key in baseline if baseline[key] != candidate[key]]
    assert len(changed) == (13 if count == 40 else 2)
    if count == 40:
        assert all(baseline[key] == candidate[key] for key in baseline if key.endswith('Program'))
        # All non-array ABI constants, flags and row shapes are unchanged.
        strip = lambda text: re.sub(r'inline constexpr std::array<std::uint32_t, \d+> (\w+)\{\{.*?\}\};', r'ARRAY \1', text, flags=re.S)
        assert strip((C / name).read_text()) == strip((B / 'stock' / name).read_text()) == strip((B / 'candidate' / name).read_text())
    headers[name] = {'arrays': count, 'changed': changed, 'stock': identity(B / 'stock' / name), 'candidate': identity(B / 'candidate' / name)}
for name in ['retained_owner_layout.h', 'retained_owner_layout.slang']:
    assert (C / name).read_bytes() == (B / 'stock' / name).read_bytes() == (B / 'candidate' / name).read_bytes()
target = W / 'products-acceptance.json'
assert not target.exists()
target.write_text(json.dumps({'passed': True, 'original_builder_owners': owners, 'modules': rows, 'headers': headers,
    'all_seven_program_tables_and_owner_layouts_unchanged': True,
    'scope': 'All35 complete modules through unchanged original builder policy; stock canonical bytes, candidate scope and supporting ABI/barrier inventories. Native/FMA and pair-camera products byte-exact. No provider, controls, benefit or scientific qualification verdict.'}, indent=2) + '\n')
print('Accepted all35 products; fifteen changed modules, complete tables/ABI and unchanged native/FMA/pair products.')
