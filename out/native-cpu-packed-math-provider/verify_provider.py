"""Bind the composed provider to preserved inputs without running kernels."""
from pathlib import Path
import bisect
import hashlib
import json
import re
import struct
import subprocess
import sys
import types
import zipfile

R = Path.cwd()
W = Path(__file__).resolve().parent
P = W / 'provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
assert W == R / 'out/native-cpu-packed-math-provider'
assert hashlib.sha256((R / 'out/native-cpu-v4-raw-provider/verify_provider.py').read_bytes()).hexdigest() == 'aa7f5c29b62d8c7f176c434624ddbfa690ec1bdd295d132a34ba7cabc61f0519'
assert hashlib.sha256((R / 'out/native-cpu-v4-provider/verify_provider.py').read_bytes()).hexdigest() == 'e60ea2063152bc100a5e738bf73fa335d8f44380385f73fd67f7d24a54b3335b'
from component_gate import component, contract_inputs, terminal
raw_component = component()
terminal(P / 'build/owner.json')
sys.path.insert(0, str(R / 'out/native-cpu-field-compatibility'))
from elf_proof import Elf, associate, pin


def doc(path):
    return json.loads(path.read_text())


def checked(path, row):
    assert pin(path) == {key: row[key] for key in ('bytes', 'sha256')}, path


def data_pin(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def named(elf, name):
    rows = [row for row in elf.symbols if row['name'] == name and row['section']]
    assert len(rows) == 1, (name, len(rows))
    return rows[0]


build, contract = doc(P / 'provider-build-result.json'), doc(P / 'source-contract.json')
old = doc(R / 'out/native-cpu-v4-raw-provider/provider/source-contract.json')
contract_inputs(contract)
assert contract['guard_shared_byte_authority'] == old['guard_shared_byte_authority']
assert contract['adapter_input_objects'] == raw_component['parts']
assert contract['raw_component_manifest'] == str(raw_component['manifest'].relative_to(R))
assert contract['math_component_manifest'] == str(raw_component['math_manifest'].relative_to(R))
assert build['passed'] and build['native_objects'] == contract['native_objects']
objects = contract['native_objects']
assert len(objects) == 8 and objects[:7] == old['native_objects'][:7]
assert contract['basis_source_identity'] == old['full_candidate_source_identity']
assert contract['new_native_compilations'] == 0 and contract['reused_raw_objects'] == 3
assert set(contract['modified_files']) == {
    'src/gallium/targets/lavapipe/meson.build',
    'src/gallium/drivers/llvmpipe/lp_screen.c',
    'src/gallium/frontends/lavapipe/lvp_device.c'}
for name, row in contract['modified_files'].items():
    checked(S / name, row['candidate'])
checked(P / 'provider.patch', contract['patch'])
unchanged_sources = {
    'src/gallium/drivers/llvmpipe/lp_native_cpu_bridge.h': dict(bytes=6421, sha256='ffa014c8d06e66c33935477fa11a13b775b9cc479860a3be5d6943adc7c3568f'),
    'src/gallium/drivers/llvmpipe/lp_state_cs.c': dict(bytes=109640, sha256='34b1f093c09543c7f1c48873081d014849088b3839316dc6b3671d2b6add9cd1')}
for name, expected in unchanged_sources.items():
    checked(S / name, expected)
guard_source = (S / 'src/gallium/drivers/llvmpipe/lp_native_cpu_bridge.h').read_text()
assert 'return util_get_cpu_caps()->has_avx2 && __builtin_cpu_supports("x86-64-v3");' in guard_source
library = R / build['library']
checked(library, build['library_identity'])
for row in objects + contract['adapter_input_objects']:
    checked(R / row['path'], row)

manifest_path = R / contract['component_manifest']
component = doc(manifest_path)
assert component['finite_component_passed'] and not component['provider_adopted']
archive = R / component['archive']['path']
checked(archive, component['archive'])
with zipfile.ZipFile(archive) as zipped:
    def member(name):
        path = 'out/native-cpu-v4-arithmetic/' + name
        data = zipped.read(path)
        assert data_pin(data) == component['members'][path]
        return json.loads(data)
    component_proof = member('code-proof.json')
    component_compile = member('compile-result.json')
assert component_proof['passed']
parts = contract['adapter_input_objects']
assert len(parts) == 7
by_name = {Path(row['path']).name: row for row in parts}
assert set(by_name) == {'rows_adapter.o', 'transport.o', 'endpoint.o', 'dopri_phase.o', 'optional.o', 'dispatch.o', 'v3.o'}
original_adapter = raw_component['adapter']
assert by_name['rows_adapter.o'] == original_adapter
# Only v3 and the baseline dispatcher retain AO object identity. The optional
# object is bound to its separate zero-residual36-site proof and finite bundle.
for label in ('dispatch', 'v3'):
    row = by_name[label + '.o']
    assert {key: row[key] for key in ('bytes', 'sha256')} == component_compile['objects'][label]
    assert component_proof['ISA_partitions'][label]['object'] == component_compile['objects'][label]
    assert component_proof['ISA_partitions'][label]['EVEX_instructions'] == 0
math_code = raw_component['math_code']
assert math_code['passed'] and math_code['no_FMA']
assert math_code['wide_product_sites'] == math_code['retention_multiply_sites'] == len(math_code['retention_dependencies']) == 36
assert {key: by_name['optional.o'][key] for key in ('bytes', 'sha256')} == math_code['optional_object']
assert math_code['original_v3_and_baseline_dispatcher_byte_preserved']

linked = Elf(library)
merged = R / objects[-1]['path']
dispatcher_path = R / by_name['dispatch.o']['path']
dispatcher = Elf(dispatcher_path)
dispatch_names = {row['name'] for row in dispatcher.symbols if (row['info'] & 15) == 2 and row['section'] and row['size']}
assert len(dispatch_names) == 11
omitted = []


class NativeElf(Elf):
    """Retain the established one-byte noop omission; dispatch is joined below."""
    def __init__(self, path):
        super().__init__(path)
        if self.kind != 1:
            return
        kept = []
        for row in self.symbols:
            if row['name'].startswith('_GLOBAL__sub_I_') and row['section'] and row['size']:
                assert self.function_bytes(row) == b'\xc3'
                assert not any(other['name'] == row['name'] and other['section'] for other in linked.symbols)
                omitted.append(dict(object=str(path), symbol=row['name'], bytes=1))
            elif Path(path).resolve() == merged.resolve() and row['name'] in dispatch_names:
                continue
            else:
                kept.append(row)
        self.symbols = kept


native_associate = types.FunctionType(associate.__code__, dict(associate.__globals__, Elf=NativeElf))
association = native_associate(library, [dict(row, path=str(R / row['path'])) for row in objects])
merge_association = associate(merged, [dict(row, path=str(R / row['path'])) for row in parts])
part_association = associate(library, [dict(row, path=str(R / row['path'])) for row in parts if Path(row['path']).name != 'dispatch.o'])

# Dispatcher joins allow only the specifically verified GOTPCRELX relaxation.
# Unlike an arbitrary opcode mask, both forms must resolve the real CPU object.
cpu = named(linked, '__cpu_features2')
dispatch_rows = []
for symbol in dispatcher.symbols:
    if symbol['name'] not in dispatch_names or not symbol['section'] or not symbol['size']:
        continue
    final = named(linked, symbol['name'])
    assert final['size'] == symbol['size']
    before = bytearray(dispatcher.function_bytes(symbol))
    after = bytearray(linked.function_bytes(final))
    assert before[:7] == bytes.fromhex('f30f1efa488b05')
    assert after[:7] == bytes.fromhex('f30f1efa488d05')
    assert final['value'] + 11 + struct.unpack_from('<i', after, 7)[0] == cpu['value']
    after[5] = before[5]
    expected = [(7, 42, -4, '__cpu_features2')]
    targets = []
    name = symbol['name']
    if name == 'sirius_native_v4_supported':
        assert before[11:] == bytes.fromhex('0fb64008c0e80283e001c3')
    else:
        assert len(before) == 37 and before[11:18] == bytes.fromhex('f6400804740fe9')
        assert before[22:32] == bytes.fromhex('662e0f1f840000000000') and before[32] == 0xe9
        for offset, version in ((18, 'v4'), (33, 'v3')):
            target = name.replace('sirius_native_', 'sirius_native_' + version + '_')
            expected.append((offset, 4, -4, target))
            address = final['value'] + offset + 4 + struct.unpack_from('<i', after, offset)[0]
            assert address == named(linked, target)['value']
            targets.append(dict(target=target, address=address))
    relocations = [(offset - symbol['value'], kind, addend, target['name'])
        for section, offset, kind, addend, target in dispatcher.relocations
        if section == symbol['section'] and 0 <= offset - symbol['value'] < symbol['size']]
    assert sorted(relocations) == sorted(expected)
    for offset, _kind, _addend, _name in relocations:
        before[offset:offset + 4] = after[offset:offset + 4] = b'\0' * 4
    assert before == after
    dispatch_rows.append(dict(symbol=name, bytes=symbol['size'], linked_address=final['value'],
        checked_GOTPCRELX_mov_to_lea=True, actual_branch_targets=targets))
assert len(dispatch_rows) == 11
association['objects'][-1]['functions'] += dispatch_rows
part_association['objects'].append(dict(object=by_name['dispatch.o'],
    all_retained_function_opcodes_match=True, functions=dispatch_rows,
    checked_GOTPCRELX_mov_to_lea=True))

# Keep the eight-row consumer schema while disclosing the exact relaxation.
rows = []
for declared, row in zip(objects, association['objects']):
    rows.append(dict(object=declared, all_retained_function_opcodes_match=True,
        functions=row['functions'], checked_dispatcher_relaxation=declared == objects[-1]))
optional_functions = next(row['functions'] for row in part_association['objects']
                          if Path(row['object']['path']).name == 'optional.o')
raw_v4_objects = []
for expected in raw_component['raw_objects']:
    found, = [row for row in part_association['objects'] if Path(row['object']['path']).resolve() == (R / expected['path']).resolve()]
    raw_v4_objects.append(dict(object=expected, all_retained_function_opcodes_match=True, functions=found['functions']))
assert [Path(row['object']['path']).stem for row in raw_v4_objects] == ['transport', 'endpoint', 'dopri_phase']
v4_ranges = {(row['linked_address'], row['linked_address'] + row['bytes']) for row in optional_functions}
v4_ranges |= {(row['linked_address'], row['linked_address'] + row['bytes']) for obj in raw_v4_objects for row in obj['functions']}
ranges = sorted((row['linked_address'], row['linked_address'] + row['bytes'],
                 (row['linked_address'], row['linked_address'] + row['bytes']) in v4_ranges)
                for obj in rows for row in obj['functions'])
starts = [start for start, _end, _v4 in ranges]
assembly = subprocess.check_output(['objdump', '-d', '--insn-width=16',
    '--start-address=' + str(ranges[0][0]), '--stop-address=' + str(max(end for _start, end, _v4 in ranges)), str(library)], text=True)
EVEX = 0
for line in assembly.splitlines():
    fields = line.split('\t')
    if len(fields) < 3 or not re.match(r'^\s*[0-9a-f]+:$', fields[0]):
        continue
    address = int(fields[0].strip()[:-1], 16)
    index = bisect.bisect_right(starts, address) - 1
    if index < 0 or address >= ranges[index][1]:
        continue
    encoded, instruction = fields[1].split(), fields[2].strip()
    assert encoded and instruction
    mnemonic = instruction.split()[0]
    assert not mnemonic.startswith(('vfm', 'vfnm', 'f')), (address, mnemonic)
    if encoded[0] == '62':
        assert ranges[index][2], (address, instruction)
        EVEX += 1
assert EVEX
blocks = {}
for block in re.split(r'\n(?=[0-9a-f]+ <)', assembly):
    match = re.match(r'([0-9a-f]+) <([^>]+)>:', block)
    if match:
        blocks[int(match[1], 16)] = block
raw_targets = {}
for prefix, target, count in (
    ('RTAddStore_0___', 'sirius_native_tagged_add8', 3),
    ('RTMultiplyStore_0___', 'sirius_native_tagged_multiply8', 3),
    ('RTSubtractStore_0___', 'sirius_native_tagged_subtract8', 3),
    ('RTDivideStore_0___', 'sirius_native_tagged_divide8', 3),
    ('RTSqrtStore_0___', 'sirius_native_tagged_sqrt8', 2)):
    found = [row for obj in rows[4:7] for row in obj['functions'] if row['symbol'].startswith(prefix)]
    assert len(found) == count
    assert all('<' + target + '>' in blocks[row['linked_address']] for row in found)
    raw_targets[target] = [dict(symbol=row['symbol'], address=row['linked_address']) for row in found]
refinement = {}
for obj in rows[4:7]:
    refinement[obj['object']['path']] = {}
    for target in ('sirius_native_dense_divide8', 'sirius_native_dense_sqrt8'):
        found = [row for row in obj['functions'] if '<' + target + '>' in blocks[row['linked_address']]]
        assert found
        refinement[obj['object']['path']][target] = [dict(symbol=row['symbol'], address=row['linked_address']) for row in found]


def disassemble(name):
    return subprocess.check_output(['objdump', '-d', '--disassemble=' + name, str(library)], text=True)


callback = disassemble('cs_exec_fn')
packet_mode = disassemble('lp_native_cpu_packet_mode')
assert '<lp_native_cpu_packet_mode>' in callback and '<util_fpstate_get>' in callback
assert callback.count('<util_fpstate_set>') == 2
assert '$0xffff003f' in callback and '$0x1f80' in callback
assert re.search(r'\bcallq?\s+\*', callback), 'Ordinary indirect JIT fallback missing'
for slug in ('transport', 'endpoint', 'dopri_phase'):
    assert '<sirius_native_' + slug + '_admit_packet>' in packet_mode
    assert '<sirius_native_' + slug + '_rows_try_packet>' in packet_mode
    wrapper = disassemble('sirius_native_' + slug + '_rows_try_packet')
    for target in ('sirius_native_' + slug + '_admit_packet', 'sirius_native_' + slug + '_try_packet',
                   'sirius_simd_' + slug + '_packet', 'sirius_native_simd_supported'):
        assert '<' + target + '>' in wrapper
    for target in ('sirius_native_v4_supported', 'sirius_simd_v4_' + slug + '_packet'):
        assert '<' + target + '>' in wrapper
for name, instruction in (('util_fpstate_get', 'stmxcsr'), ('util_fpstate_set', 'ldmxcsr')):
    assert instruction in disassemble(name)
native_guard = disassemble('sirius_native_simd_supported')
assert '<__cpu_features2>' in native_guard and re.search(r'shr\s+\$1,%al', native_guard)
assert 'testb' in native_guard and '<call_once>' in native_guard
batch = doc(R / 'out/native-cpu-slang-bridge/batch-generation.json')
modules = {}
for stage in ('Transport', 'Endpoint', 'DopriPhase'):
    value = batch['modules'][stage]['SPIRV']
    checked(R / value['module'], value)
    assert linked.data.count((R / value['module']).read_bytes()) == 1
    modules[stage] = {key: value[key] for key in ('module', 'bytes', 'sha256')}
notes = subprocess.check_output(['readelf', '-n', str(library)], text=True)
needed = re.findall(r'x86 ISA needed:\s*([^\n]+)', notes)
assert not needed or needed == ['x86-64-baseline'], needed
resolution = subprocess.check_output(['ldd', '-r', str(library)], text=True, stderr=subprocess.STDOUT)
assert 'undefined symbol' not in resolution and 'not found' not in resolution
# Join every declared application call in new raw objects and the selector.
# These are exact PLT32 fields, resolving to defined native functions here.
def direct_joins(path, functions, keep):
    obj = Elf(path)
    joins = []
    for section, location, kind, addend, target in obj.relocations:
        if not keep(target['name']):
            continue
        assert kind == 4 and addend == -4, (path, target['name'], kind, addend)
        owners = [s for s in obj.symbols if s['section'] == section and (s['info'] & 15) == 2
                  and s['size'] and s['value'] <= location < s['value'] + s['size']]
        assert owners
        addresses = set()
        for owner in owners:
            row, = [f for f in functions if f['symbol'] == owner['name']]
            offset = location - owner['value']
            final, = [s for s in linked.symbols if s['name'] == owner['name'] and s['section']
                      and s['value'] == row['linked_address'] and s['size'] == row['bytes']]
            code = linked.function_bytes(final)
            assert 0 < offset and offset + 4 <= len(code) and code[offset - 1] in (0xe8, 0xe9)
            actual = final['value'] + offset + 4 + struct.unpack_from('<i', code, offset)[0]
            expected = named(linked, target['name'])['value']
            assert actual == expected
            addresses.add(actual)
        assert len(addresses) == 1
        joins.append(dict(caller=owners[0]['name'], target=target['name'], linked_address=addresses.pop(),
                          relocation_type=kind, addend=addend))
    assert joins
    return joins

raw_direct_helper_joins = {}
for item in raw_v4_objects:
    path = R / item['object']['path']
    joins = direct_joins(path, item['functions'], lambda n: n.startswith('sirius_native_v4_'))
    assert {row['target'] for row in joins} == set(next(row['undefined'] for row in raw_component['code']['objects']
        if row['candidate'] == item['object'])) - {'memset', 'sirius_simd_fmaxf'}
    raw_direct_helper_joins[item['object']['path']] = joins
adapter_rows, = [row for row in part_association['objects']
                 if Path(row['object']['path']).resolve() == (R / by_name['rows_adapter.o']['path']).resolve()]
selector_joins = direct_joins(R / by_name['rows_adapter.o']['path'], adapter_rows['functions'],
    lambda n: n.startswith('sirius_native_') or n.startswith('sirius_simd_') and n.endswith('_packet'))
for slug in ('transport', 'endpoint', 'dopri_phase'):
    found = {row['target'] for row in selector_joins if row['caller'] == 'sirius_native_' + slug + '_rows_try_packet'}
    assert found == {'sirius_native_simd_supported', 'sirius_native_v4_supported',
        'sirius_native_' + slug + '_try_packet', 'sirius_native_' + slug + '_admit_packet',
        'sirius_simd_' + slug + '_packet', 'sirius_simd_v4_' + slug + '_packet'}
# Count the exact retained functions in these input objects, rather than
# carrying forward the predecessor's function count across changed code.
verified_omissions = {(str(Path(row['object']).resolve()), row['symbol']) for row in omitted}
expected_functions = 0
for item in objects:
    path = R / item['path']
    obj = Elf(path)
    for symbol in obj.symbols:
        if (symbol['info'] & 15) != 2 or not symbol['section'] or not symbol['size']:
            continue
        if (str(path.resolve()), symbol['name']) not in verified_omissions:
            expected_functions += 1
assert sum(len(row['functions']) for row in rows) == expected_functions

result = dict(passed=True, library=build['library'], library_identity=pin(library), objects=rows,
    first7_native_objects_unchanged=True, original_v3_component_proof_reused=True,
    raw_v4_objects=raw_v4_objects, raw_v4_direct_helper_joins=raw_direct_helper_joins,
    baseline_selector_function_joins=selector_joins,
    raw_component_manifest=dict(path=contract['raw_component_manifest'], **pin(raw_component['manifest'])),
    guard_shared_byte_authority=contract['guard_shared_byte_authority'],
    component_manifest=dict(path=contract['component_manifest'], **pin(manifest_path)),
    component_code_proof=component['members']['out/native-cpu-v4-arithmetic/code-proof.json'],
    math_component_manifest=dict(path=contract['math_component_manifest'], **pin(raw_component['math_manifest'])),
    optional_math_code_proof=raw_component['math_code_proof'],
    checked_optional_math_splitter_retention_dependencies=36,
    expected_retained_input_functions=expected_functions,
    unchanged_guard_and_callback_sources=unchanged_sources,
    seven_input_merged_code=merge_association, seven_input_final_code=part_association,
    dispatcher_function_joins=dispatch_rows, omitted_exact_one_byte_RET_initializers=omitted,
    callback_shared_mode_reaches_rows_pure_SIMD_and_scalar_all3=True,
    callback_saved_and_restored_MXCSR_instructions_present=True,
    original_v3_CPU_guard_and_ordinary_fallback_retained=True,
    native_FMA_mnemonics=0, optional_EVEX_instructions=EVEX,
    EVEX_confined_to_separately_guarded_v4_functions=True,
    all3_complete_SPIRV_spans_present_once=True, complete_modules=modules,
    dense_C_ABI_actual_targets_all3=raw_targets, complete_refinement_CPP_targets_all3=refinement,
    linked_runtime_symbols_resolve=True, runtime_resolution=resolution, linked_ELF_notes=notes,
    verifier=pin(Path(__file__)),
    scope='Static source/ELF seven-input AR composition only. First7 AP objects remain exact; newraw3 leaf joins are separate from the eight-row schema. '
          'dispatcher code permits only the checked GOTPCRELX mov-to-lea relaxation to the actual CPU symbol. '
          'All other declared relocation operands retain the existing four-byte masking bounds; matching '
          'function bytes alone is not a universal proof of relocation or external-library semantics. '
          'The original v3/dispatcher ISA evidence is reused; new optional math is bound to its sealed36-site proof. '
          'all final native EVEX instructions are confined to the sealed guarded raw-v4 or math-v4 ranges. ' 
          'No provider loading/execution, native runtime admission, numerical correctness, older-hardware '
          'fallback execution, benefit/frame acceptance or scientific/release qualification is established.')
with (P / 'linked-native-proof.json').open('x') as stream:
    stream.write(json.dumps(result, indent=2) + '\n')
print(json.dumps(dict(passed=True, native_objects=len(rows),
    linked_functions=sum(len(row['functions']) for row in rows), dispatcher_joins=len(dispatch_rows),
    complete_SPIRV_modules=len(modules), runtime_symbols_resolve=True)))
