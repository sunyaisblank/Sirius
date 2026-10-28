"""AU verifier derived from pinned AT verifier 5ba47647; no workload execution."""
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
assert W == R / 'out/native-cpu-shared-value-provider'
assert hashlib.sha256((R / 'out/native-cpu-dispatch-token-provider/verify_provider.py').read_bytes()).hexdigest() == '5ba476473f420d9dfdbcc056593336142b98d608f0995afe91db3617ec473009'
assert hashlib.sha256((R / 'out/native-cpu-inline-private-provider/verify_provider.py').read_bytes()).hexdigest() == '920a1a71e0ddd3834e45abd0cea264909d8f3c265a318b58985348e47860aa68'
assert hashlib.sha256((R / 'out/native-cpu-packed-math-provider/verify_provider.py').read_bytes()).hexdigest() == 'b70ca03e5450cb876baddfaa5e6e611375a86e650d97d254aff70cf2a5d651a6'
assert hashlib.sha256((R / 'out/native-cpu-v4-raw-provider/verify_provider.py').read_bytes()).hexdigest() == 'aa7f5c29b62d8c7f176c434624ddbfa690ec1bdd295d132a34ba7cabc61f0519'
assert hashlib.sha256((R / 'out/native-cpu-v4-provider/verify_provider.py').read_bytes()).hexdigest() == 'e60ea2063152bc100a5e738bf73fa335d8f44380385f73fd67f7d24a54b3335b'
from component_gate import component, contract_inputs, terminal, HEADERS, CORES
raw_component = component()
# Historical finite evidence justifies only the unchanged public/legacy objects.
# It never establishes numerical qualification of the new shared interface.
import importlib.util
legacy_gate_path = R / 'out/native-cpu-dispatch-token-provider/component_gate.py'
assert hashlib.sha256(legacy_gate_path.read_bytes()).hexdigest() == '73f6a7b31f539b7d70df30f57f0824202c8cb9b4fd869147ff6e455f25e41875'
spec = importlib.util.spec_from_file_location('AU_legacy_AT_gate', legacy_gate_path)
legacy_gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy_gate)
legacy = legacy_gate.component()

terminal(P / 'build/owner.json')
sys.path.insert(0, str(R / 'out/native-cpu-field-compatibility'))
from elf_proof import Elf as OriginalElf, associate, pin


class Elf(OriginalElf):
    """Cache immutable section/function spans during exact repeated joins."""
    def __init__(self, path):
        self._section_bytes = {}
        self._function_bytes = {}
        super().__init__(path)

    def bytes(self, index):
        if index not in self._section_bytes:
            self._section_bytes[index] = super().bytes(index)
        return self._section_bytes[index]

    def function_bytes(self, row):
        key = (row['section'], row['value'], row['size'])
        if key not in self._function_bytes:
            self._function_bytes[key] = super().function_bytes(row)
        return self._function_bytes[key]


# The inherited associator must use the same cached reader; every comparison
# and relocation bound is otherwise the exact sealed implementation.
associate = types.FunctionType(associate.__code__, dict(associate.__globals__, Elf=Elf))


def defined(elf, name):
    if not hasattr(elf, '_defined_by_name'):
        elf._defined_by_name = {}
        for row in elf.symbols:
            if row['section']:
                elf._defined_by_name.setdefault(row['name'], []).append(row)
    return elf._defined_by_name.get(name, [])


def doc(path):
    return json.loads(path.read_text())


def checked(path, row):
    assert pin(path) == {key: row[key] for key in ('bytes', 'sha256')}, path


def data_pin(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def named(elf, name):
    rows = defined(elf, name)
    assert len(rows) == 1, (name, len(rows))
    return rows[0]


build, contract = doc(P / 'provider-build-result.json'), doc(P / 'source-contract.json')
old = doc(R / 'out/native-cpu-dispatch-token-provider/provider/source-contract.json')
contract_inputs(contract)
assert contract['guard_shared_byte_authority'] == old['guard_shared_byte_authority'] == legacy['guard_shared_byte_authority']
assert contract['adapter_input_objects'] == raw_component['parts']
assert contract['component_interface'] == raw_component['component_interface']
assert contract['shared_value_interface'] == raw_component['shared_value_interface']
assert not contract['component_interface']['numerical_qualification_transferred']
assert contract['component_interface']['genuine_job_numerical_gate_required']
assert build['passed'] and build['native_objects'] == contract['native_objects']
objects = contract['native_objects']
assert len(objects) == 8 and objects[:7] == raw_component['native_first7'] == old['native_objects'][:7]
assert objects[-1] == raw_component['merged']
assert contract['basis_source_identity'] == old['full_candidate_source_identity']
assert contract['new_native_compilations'] == 0
assert contract['provider_identity'] == 'Sirius-native-cpu-shared-value-v45'
assert contract['provider_UUID'] == 'SIR-CPU-25.2.8AU'
assert set(contract['modified_files']) == {
    'src/gallium/targets/lavapipe/meson.build',
    'src/gallium/drivers/llvmpipe/lp_screen.c',
    'src/gallium/frontends/lavapipe/lvp_device.c'}
for name, row in contract['modified_files'].items():
    checked(S / name, row['candidate'])
checked(P / 'provider.patch', contract['patch'])
unchanged_sources = {
    'src/gallium/drivers/llvmpipe/lp_state_cs.c': dict(bytes=109640, sha256='34b1f093c09543c7f1c48873081d014849088b3839316dc6b3671d2b6add9cd1')}
for name, expected in unchanged_sources.items():
    checked(S / name, expected)
candidate_headers = {name: (S / name).read_text() for name in HEADERS}
for name, text in candidate_headers.items():
    expected = old['modified_files'][name]['candidate']
    assert data_pin(text.encode()) == expected == raw_component['unchanged_bridge_headers'][name]
    assert contract['source_inputs'][str((S / name).relative_to(R))] == expected
assert contract['observer_layout_assertion_identity'] == old['observer_layout_assertion_identity'] == data_pin(legacy_gate.LAYOUT_ASSERTIONS.encode())
restored_headers = legacy_gate.header_patch(legacy_gate.remove_layout_assertions(candidate_headers), legacy['mesa_patch'], reverse=True)
for name, text in restored_headers.items():
    assert data_pin(text.encode()) == legacy['header_basis'][name]
assert legacy_gate.add_layout_assertions(legacy_gate.header_patch(restored_headers, legacy['mesa_patch'])) == candidate_headers
assert contract['admission_header_patch_identity'] == old['admission_header_patch_identity'] == data_pin(legacy['mesa_patch'].encode())
assert contract['admission_component_manifest'] == old['admission_component_manifest'] == str(legacy['manifest'].relative_to(R))
guard_source = candidate_headers[HEADERS[0]]
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
assert len(parts) == 11
by_name = {Path(row['path']).name: row for row in parts}
assert len(by_name) == 11 and set(by_name) == {
    'adapter-only.o', 'transport.o', 'endpoint.o', 'dopri_phase.o', 'optional.o', 'dispatch.o', 'v3.o',
    'shared-dp.o', 'shared-transport.o', 'shared-endpoint.cleanup.o', 'shared-math.o'}
assert by_name['adapter-only.o'] == raw_component['parts'][0]
assert parts[1:7] == legacy['parts'][1:]
for label in ('dispatch', 'v3'):
    row = by_name[label + '.o']
    assert {key: row[key] for key in ('bytes', 'sha256')} == component_compile['objects'][label]
    assert component_proof['ISA_partitions'][label]['object'] == component_compile['objects'][label]
    assert component_proof['ISA_partitions'][label]['EVEX_instructions'] == 0
math_code = legacy['math_code']
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

# Exact public function bytes and all call relocations stay those of AT.
# The separately named private callbacks alone can reach the new scoped raws.
public_wrapper_byte_inverse = []
old_adapter = Elf(R / legacy['adapter']['path'])
new_adapter = Elf(R / by_name['adapter-only.o']['path'])
for stage in ('transport', 'endpoint', 'dopri_phase'):
    name = f'sirius_native_{stage}_rows_try_packet'
    before, after = named(old_adapter, name), named(new_adapter, name)
    assert before['size'] == after['size'] == 295
    codes = [bytearray(obj.function_bytes(row)) for obj, row in ((old_adapter, before), (new_adapter, after))]
    relocs = []
    for obj, row, code in ((old_adapter, before, codes[0]), (new_adapter, after, codes[1])):
        operands = []
        for section, location, kind, addend, target in obj.relocations:
            offset = location - row['value']
            if section == row['section'] and 0 <= offset < row['size']:
                assert kind == 4 and addend == -4 and 0 < offset and offset + 4 <= len(code)
                assert code[offset - 1] in (0xe8, 0xe9)
                operands.append((offset, kind, addend, target['name']))
                code[offset:offset + 4] = b'\0' * 4
        relocs.append(operands)
    assert len(relocs[0]) == 6 and relocs[0] == relocs[1] and codes[0] == codes[1]
    public_wrapper_byte_inverse.append(dict(symbol=name, bytes=295,
        nonrelocation_bytes_identical=True, exact_call_relocations=relocs[0]))

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
def associated(expected):
    found, = [row for row in part_association['objects']
              if Path(row['object']['path']).resolve() == (R / expected['path']).resolve()]
    return dict(object=expected, all_retained_function_opcodes_match=True, functions=found['functions'])

raw_v4_objects = [associated(row) for row in raw_component['raw_v4_objects']]
old_public_raw_v4_objects = [associated(row) for row in legacy['raw_objects']]
shared_value_math_object = associated(raw_component['shared_value_math_object'])
math_object = Elf(R / shared_value_math_object['object']['path'])
assert {row['name'] for row in math_object.symbols if row['name'].startswith('sirius_shared_value_') and row['section']} == CORES
for name in CORES:
    symbol = named(math_object, name)
    assert (symbol['info'] & 15) == 2 and symbol['info'] >> 4 == 1 and symbol['other'] & 3 == 2 and symbol['size']
    final = named(linked, name)
    assert final['other'] & 3 == 2 or final['info'] >> 4 == 0
    assert any(row['symbol'] == name and row['linked_address'] == final['value'] and row['bytes'] == final['size']
               for row in shared_value_math_object['functions'])
for item, scope in zip(raw_v4_objects, ('transport', 'endpoint', 'dp')):
    assert any(row['symbol'] == 'sirius_scope_' + scope + '_packet' for row in item['functions'])
new_value_ranges = {(row['linked_address'], row['linked_address'] + row['bytes'])
                    for obj in raw_v4_objects + [shared_value_math_object] for row in obj['functions']}
v4_ranges = {(row['linked_address'], row['linked_address'] + row['bytes']) for row in optional_functions}
v4_ranges |= new_value_ranges
v4_ranges |= {(row['linked_address'], row['linked_address'] + row['bytes']) for obj in old_public_raw_v4_objects for row in obj['functions']}
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
    assert mnemonic not in ('(bad)', '.byte'), (address, instruction)
    if (ranges[index][0], ranges[index][1]) in new_value_ranges:
        assert not mnemonic.startswith(('ldmxcsr', 'stmxcsr', 'xsave', 'xrstor', 'fxsave', 'fxrstor')), (address, instruction)
        assert not (mnemonic in ('call', 'callq') and instruction.split(None, 1)[1].startswith('*')), (address, instruction)
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
    assert '<sirius_native_' + slug + '_rows_admitted_packet>' in packet_mode
    wrapper = disassemble('sirius_native_' + slug + '_rows_try_packet')
    for target in ('sirius_native_' + slug + '_admit_packet', 'sirius_native_' + slug + '_try_packet',
                   'sirius_simd_' + slug + '_packet', 'sirius_native_simd_supported'):
        assert '<' + target + '>' in wrapper
    for target in ('sirius_native_v4_supported', 'sirius_simd_v4_' + slug + '_packet'):
        assert '<' + target + '>' in wrapper
    private = 'sirius_native_' + slug + '_rows_admitted_packet'
    private_wrapper = disassemble(private)
    for target in ('sirius_native_simd_supported', 'sirius_native_v4_supported',
                   'sirius_native_' + slug + '_admitted_packet',
                   'sirius_simd_' + slug + '_packet', 'sirius_scope_' + ('dp' if slug == 'dopri_phase' else slug) + '_packet'):
        assert '<' + target + '>' in private_wrapper
    for forbidden in ('memcmp', 'bcmp', 'sirius_native_' + slug + '_admit_packet',
                      'sirius_native_' + slug + '_try_packet'):
        assert '<' + forbidden + '>' not in private_wrapper
    assert named(Elf(R/raw_component['parts'][0]['path']), private)['other'] & 3 == 2
    scalar_private = 'sirius_native_' + slug + '_admitted_packet'
    assert named(Elf(R/objects[('transport','endpoint','dopri_phase').index(slug)]['path']), scalar_private)['other'] & 3 == 2
    for private_name in (private, scalar_private):
        symbol=named(linked,private_name)
        assert symbol['other'] & 3 == 2 or symbol['info'] >> 4 == 0
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
            final, = [s for s in defined(linked, owner['name'])
                      if s['value'] == row['linked_address'] and s['size'] == row['bytes']]
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

raw_direct_helper_joins, old_public_direct_helper_joins = {}, {}
for item, stage, count in zip(raw_v4_objects, ('transport', 'endpoint', 'dopri_phase'), (2701, 3959, 1581)):
    joins = direct_joins(R / item['object']['path'], item['functions'], lambda n: n.startswith('sirius_shared_value_'))
    expected = CORES - ({'sirius_shared_value_tagged_sqrt8'} if stage == 'dopri_phase' else set())
    assert len(joins) == count and {row['target'] for row in joins} == expected
    raw_direct_helper_joins[item['object']['path']] = joins
assert sum(map(len, raw_direct_helper_joins.values())) == 8241
for item in old_public_raw_v4_objects:
    joins = direct_joins(R / item['object']['path'], item['functions'], lambda n: n.startswith('sirius_native_v4_'))
    assert {row['target'] for row in joins} == set(next(row['undefined'] for row in legacy['code']['objects']
        if row['candidate'] == item['object'])) - {'memset', 'sirius_simd_fmaxf'}
    old_public_direct_helper_joins[item['object']['path']] = joins
adapter_rows, = [row for row in part_association['objects']
                 if Path(row['object']['path']).resolve() == (R / by_name['adapter-only.o']['path']).resolve()]
selector_joins = direct_joins(R / by_name['adapter-only.o']['path'], adapter_rows['functions'],
    lambda n: n.startswith('sirius_native_') or n.startswith(('sirius_simd_', 'sirius_scope_')) and n.endswith('_packet'))
for slug in ('transport', 'endpoint', 'dopri_phase'):
    found = {row['target'] for row in selector_joins if row['caller'] == 'sirius_native_' + slug + '_rows_try_packet'}
    assert found == {'sirius_native_simd_supported', 'sirius_native_v4_supported',
        'sirius_native_' + slug + '_try_packet', 'sirius_native_' + slug + '_admit_packet',
        'sirius_simd_' + slug + '_packet', 'sirius_simd_v4_' + slug + '_packet'}
    private_found = {row['target'] for row in selector_joins if row['caller'] == 'sirius_native_' + slug + '_rows_admitted_packet'}
    assert private_found == {'sirius_native_simd_supported', 'sirius_native_v4_supported',
        'sirius_native_' + slug + '_admitted_packet', 'sirius_simd_' + slug + '_packet',
        'sirius_scope_' + ('dp' if slug == 'dopri_phase' else slug) + '_packet'}
shared_value_interface = dict(raw_component['shared_value_interface'], no_FP_observers=True)

# Check the concrete positive-byte branch, not just the existence of private
# symbol names. The job layout is the unchanged pinned lp_state_cs.c schema;
# the real pointer/stage/views/positive flag are independently read during C8.
mode_instructions=[]
for line in packet_mode.splitlines():
    match=re.match(r'^\s*([0-9a-f]+):\s+(?:[0-9a-f]{2}\s+)+\s*([a-z][a-z0-9]*)\s*(.*)$',line)
    if match:mode_instructions.append((int(match[1],16),match[2],match[3]))
assert mode_instructions
addresses={address:index for index,(address,_name,_args) in enumerate(mode_instructions)}

def successors(index):
    _address,mnemonic,args=mode_instructions[index]
    next_index=index+1
    if mnemonic.startswith('ret') or mnemonic in ('ud2','hlt'):return []
    if mnemonic in ('call','callq') and '<__stack_chk_fail@plt>' in args:return []
    if mnemonic.startswith('j'):
        match=re.match(r'([0-9a-f]+)\s',args)
        assert match,('Unknown packet-mode branch',args)
        target=addresses[int(match[1],16)]
        return [target] if mnemonic in ('jmp','jmpq') else [target,next_index]
    assert next_index<len(mode_instructions)
    return [next_index]

def reachable(start, stop=None):
    found=set();pending=[start]
    while pending:
        index=pending.pop()
        if index==stop or index in found:continue
        found.add(index);pending.extend(successors(index))
    return found

private_mode_calls={index for index,(_address,mnemonic,args) in enumerate(mode_instructions)
    if mnemonic in ('call','callq','jmp','jmpq') and
    any('<sirius_native_'+slug+'_rows_admitted_packet>' in args for slug in ('transport','endpoint','dopri_phase'))}
assert len(private_mode_calls)==3
# These exact data reads are compile-time asserted in the derived private
# header, and are tied below to the unchanged callback's callee-saved job.
layout=dict(job_bytes=112, grid_size_offset=0, native_preflight_passed_offset=63,
    stage_offset=65,current_offset=72, table_view_offset=24,buffer_bytes=16,
    buffer_pointer_offset=0,buffer_extent_offset=8,descriptor_bytes=64)
# A small must-alias transfer binds the flag's base register to the incoming
# job pointer on every CFG path. Unsupported pointer reload/spill forms stop
# for explicit review; no memory location is guessed from a matching offset.
register_aliases={}
for full,parts in {'rax':('rax','eax','ax','al','ah'),'rbx':('rbx','ebx','bx','bl','bh'),
                  'rcx':('rcx','ecx','cx','cl','ch'),'rdx':('rdx','edx','dx','dl','dh'),
                  'rsi':('rsi','esi','si','sil'),'rdi':('rdi','edi','di','dil'),
                  'rbp':('rbp','ebp','bp','bpl'),'rsp':('rsp','esp','sp','spl')}.items():
    for name in parts:register_aliases[name]=full
for number in range(8,16):
    for suffix in ('','d','w','b'):register_aliases['r'+str(number)+suffix]='r'+str(number)
volatile={'rax','rcx','rdx','rsi','rdi','r8','r9','r10','r11'}

def job_alias_transfer(index, incoming):
    aliases=set(incoming);_address,mnemonic,args=mode_instructions[index]
    if mnemonic in ('call','callq'):aliases-=volatile
    if mnemonic.startswith(('rep','div','idiv','mul')):
        aliases-=({'rdi','rsi','rcx'} if mnemonic.startswith('rep') else {'rax','rdx'})
    self_exchange=re.fullmatch(r'%([a-z0-9]+),%\1',args) if mnemonic=='xchg' else None
    assert mnemonic != 'xchg' or (self_exchange and
        (self_exchange[1] == register_aliases.get(self_exchange[1]) or self_exchange[1] == 'ax')), 'Unsupported pointer exchange'
    assert mnemonic not in ('cpuid','rdtsc','rdtscp','lahf','enter','enterq','popa','popal','popaq','xadd','cmpxchg','cmpxchg8b','cmpxchg16b'), 'Unsupported implicit GP clobber'
    assert not mnemonic.startswith('loop'), 'Unsupported implicit counter branch'
    if mnemonic in ('leave','leaveq'):aliases -= {'rbp','rsp'}
    if mnemonic.startswith(('movs','stos','lods','cmps','scas')):aliases -= {'rax','rdi','rsi','rcx'}
    if mnemonic.startswith('imul'):aliases -= {'rax','rdx'}
    if not self_exchange and mnemonic not in ('cmp','cmpq','cmpl','cmpb','cmpw','test','testq','testl','testb','testw','push','pushq','nop','nopl','nopw'):
        destination=re.search(r'(?:^|,)%([a-z0-9]+)$',args)
        if destination and destination[1] in register_aliases:
            target=register_aliases[destination[1]]
            source=re.fullmatch(r'%([a-z0-9]+),%([a-z0-9]+)',args)
            inherited=(mnemonic in ('mov','movq') and source and source[1] in incoming
                       and source[1]==register_aliases.get(source[1])
                       and source[2]==target)
            aliases.discard(target)
            if inherited:aliases.add(target)
    return aliases

job_alias_states={0:{'rdi'}};pending=[0]
while pending:
    index=pending.pop();outgoing=job_alias_transfer(index,job_alias_states[index])
    for next_index in successors(index):
        previous=job_alias_states.get(next_index)
        merged=outgoing if previous is None else previous & outgoing
        if previous is None or merged!=previous:
            job_alias_states[next_index]=set(merged);pending.append(next_index)
positive_branches=[]
for index,(address,mnemonic,args) in enumerate(mode_instructions[:-1]):
    # lp_cs_job_info.native_preflight_passed is the byte at63 after the
    # unchanged three grid/iter/base/block arrays, three uints and three bools.
    if mnemonic!='cmpb' or not re.fullmatch(r'\$0x0,0x3f\(%[a-z0-9]+\)',args):continue
    branch=mode_instructions[index+1]
    if branch[1] not in ('je','jz','jne','jnz'):continue
    flag_base=re.search(r'\(%([a-z0-9]+)\)',args)[1]
    assert flag_base in job_alias_states[index], 'Flag base not proven equal to incoming job'
    target,fallthrough=successors(index+1)
    false,true=(target,fallthrough) if branch[1] in ('je','jz') else (fallthrough,target)
    assert not private_mode_calls & reachable(false)
    assert private_mode_calls <= reachable(true)
    assert not private_mode_calls & reachable(0,stop=index)
    positive_branches.append(dict(flag_offset=63,flag_base_register=flag_base,
        flag_base_must_alias_incoming_job=True,compare_address=address,compare=mnemonic+' '+args,
        conditional_address=branch[0],condition=branch[1],
        positive_successor=mode_instructions[true][0],negative_successor=mode_instructions[false][0],
        dominates_all_three_private_calls=True,false_branch_reaches_no_private_call=True))
assert len(positive_branches)==1, 'Review actual compiler job-capability branch before loading'
# The callback retains its genuine incoming job in RBX across the one mode
# call. The normal unwind at an existing worker-entry stop restores that
# callee-saved value. No debug locals or added runtime breakpoint is needed.
callback_instructions=[]
for line in callback.splitlines():
    match=re.match(r'^\s*([0-9a-f]+):\s+(?:[0-9a-f]{2}\s+)+\s*([a-z][a-z0-9]*)\s*(.*)$',line)
    if match:callback_instructions.append((int(match[1],16),match[2],match[3]))
mode_calls=[index for index,(_address,mnemonic,args) in enumerate(callback_instructions)
    if mnemonic in ('call','callq') and '<lp_native_cpu_packet_mode>' in args]
mode_call,=mode_calls
prefix=callback_instructions[:mode_call]
job_saves=[row for row in prefix if row[1] in ('mov','movq') and row[2]=='%rdi,%rbx']
assert len(job_saves)==1
job_writes=[row for row in prefix if re.search(r',%(?:rbx|ebx|bx|bl|bh)$',row[2])
    and row[1] not in ('cmp','cmpq','cmpl','cmpb','test','testq','testl','testb')]
assert job_writes==job_saves
assert any(mnemonic in ('mov','movq') and args=='%rbx,%rdi' for _address,mnemonic,args in prefix[-8:])
r8_writes=[row for row in prefix if re.search(r',%r8(?:d|w|b)?$',row[2])]
assert r8_writes and r8_writes[-1][1]=='xor' and r8_writes[-1][2]=='%r8d,%r8d'
after=callback_instructions[mode_call+1]
assert after[1]=='cmpb' and after[2]=='$0x0,0x3f(%rbx)'
callback_addresses={address:index for index,(address,_name,_args) in enumerate(callback_instructions)}
def callback_reachable_without(skipped):
    found=set();pending=[0]
    while pending:
        index=pending.pop()
        if index==skipped or index in found:continue
        found.add(index)
        _address,mnemonic,args=callback_instructions[index]
        if mnemonic.startswith('ret') or mnemonic in ('ud2','hlt'):continue
        if mnemonic in ('call','callq') and '<__stack_chk_fail@plt>' in args:continue
        if mnemonic.startswith('j'):
            match=re.match(r'([0-9a-f]+)\s',args)
            assert match, ('Unsupported callback branch',args)
            target=callback_addresses.get(int(match[1],16))
            if target is None:
                assert mnemonic in ('jmp','jmpq'), 'Conditional branch leaves callback'
                continue
            pending.append(target)
            if mnemonic in ('jmp','jmpq'):continue
        assert index+1<len(callback_instructions)
        pending.append(index+1)
    return found
save_index=callback_addresses[job_saves[0][0]]
assert mode_call not in callback_reachable_without(save_index), 'Incoming job save does not dominate callback call'
argument_writes=[]
for index,(_address,mnemonic,args) in enumerate(prefix):
    assert mnemonic not in ('cpuid','popa','popal','popaq'), 'Unsupported implicit callback job clobber'
    if mnemonic=='xchg':
        assert re.fullmatch(r'%([a-z0-9]+),%\1',args) and not re.search(r'%(?:rbx|ebx|bx|bl|bh)\b',args), 'Unsupported callback job exchange'
    if mnemonic in ('call','callq') or re.search(r',%(?:rdi|edi|di|dil)$',args) or mnemonic.startswith(('rep','movs','stos','lods','cmps','scas')):
        argument_writes.append(index)
argument_index=argument_writes[-1]
assert prefix[argument_index][1] in ('mov','movq') and prefix[argument_index][2]=='%rbx,%rdi', 'Final mode argument is not the saved job'
assert mode_call not in callback_reachable_without(argument_index), 'Saved-job argument does not dominate callback call'
assert not any(row[1] in ('call','callq') for row in prefix[argument_index+1:]), 'Intervening call clobbers job argument'
r8_index=callback_addresses[r8_writes[-1][0]]
assert mode_call not in callback_reachable_without(r8_index), 'False-preflight argument does not dominate callback call'
assert not any(row[1] in ('call','callq') for row in prefix[r8_index+1:]), 'Intervening call clobbers preflight argument'
callback_job=dict(frame='cs_exec_fn',register='rbx',function_address=named(linked,'cs_exec_fn')['value'],
    packet_mode_address=named(linked,'lp_native_cpu_packet_mode')['value'],incoming_job_save=job_saves[0][0],
    packet_mode_call=callback_instructions[mode_call][0], return_address=after[0],
    preflight_false_argument=r8_writes[-1][0],postcall_positive_job_read=after[0],
    no_intervening_job_register_writes=True, incoming_job_save_dominates_call=True,
    final_job_argument_definition=prefix[argument_index][0], final_job_argument_dominates_call=True,
    no_intervening_argument_call=True, false_preflight_argument_dominates_call=True,
    unwound_callee_saved_register_required=True)
mesa_worker_calls=[]
for index in sorted(private_mode_calls):
    address,mnemonic,args=mode_instructions[index]
    target=re.search(r'<([^>]+)>',args)[1]
    actual=int(re.match(r'([0-9a-f]+)',args)[1],16)
    assert actual==named(linked,target)['value']
    mesa_worker_calls.append(dict(instruction=address,mnemonic=mnemonic,target=target,actual_address=actual,
        return_address=mode_instructions[index+1][0]))
# Bind full FP status/control save, IEEE32 controls, and full restore to the
# same actual callback. Accept only the known register/operation forms.
fp_get, = [i for i, row in enumerate(callback_instructions) if row[1] in ('call','callq') and '<util_fpstate_get>' in row[2]]
fp_sets = [i for i, row in enumerate(callback_instructions) if row[1] in ('call','callq') and '<util_fpstate_set>' in row[2]]
assert len(fp_sets) == 2
fp_enter, fp_restore = fp_sets
assert [(op,args) for _pc,op,args in callback_instructions[fp_get+1:fp_enter]] == [
    ('mov','%eax,%edi'), ('mov','%eax,%r13d'), ('and','$0xffff003f,%edi'), ('or','$0x1f80,%edi')]
assert fp_enter < mode_call < fp_restore
assert callback_instructions[fp_restore-1][1:] == ('mov','%r13d,%edi')
r13_writes = [i for i, (_pc,op,args) in enumerate(callback_instructions[:fp_restore])
              if re.search(r',%(?:r13|r13d|r13w|r13b)$',args) and op not in ('cmp','test')]
assert r13_writes == [fp_get+2]
assert mode_call not in callback_reachable_without(fp_get)
assert mode_call not in callback_reachable_without(fp_enter)
# All normal hot returns after the mode call must pass the one full restore.
# The one source-bound cold edge is the existing noreturn abort, not a return.
cold = named(linked, 'cs_exec_fn.cold')
assert linked.function_bytes(cold)[0] == 0xe8 and cold['size'] == 5
cold_target = cold['value'] + 5 + struct.unpack_from('<i', linked.function_bytes(cold), 1)[0]
assert '<abort@plt>' in disassemble('cs_exec_fn.cold')
assert cold_target == int(re.search(r'call\s+([0-9a-f]+)\s+<abort@plt>', disassemble('cs_exec_fn.cold'))[1], 16)
seen, pending = set(), [mode_call+1]
while pending:
    i = pending.pop()
    if i == fp_restore or i in seen: continue
    seen.add(i)
    _pc, op, args = callback_instructions[i]
    assert not op.startswith('ret'), 'Normal callback return bypasses full FP restore'
    if op in ('call','callq') and '<__stack_chk_fail@plt>' in args: continue
    if op.startswith('j'):
        target = int(re.match(r'([0-9a-f]+)\s',args)[1],16)
        if target in callback_addresses: pending.append(callback_addresses[target])
        else: assert target == cold['value']
        if op in ('jmp','jmpq'): continue
    assert i+1 < len(callback_instructions)
    pending.append(i+1)
fp_scope = dict(saved_full_mxcsr=True, traps_masked=True, round_nearest=True,
    daz_ftz_off=True, restored_full_mxcsr=True, actual_callback_code_join=True,
    get_address=callback_instructions[fp_get][0], enter_address=callback_instructions[fp_enter][0],
    restore_address=callback_instructions[fp_restore][0], saved_register='r13d',
    source=unchanged_sources, same_positive_job_branch_required=True)
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
    unchanged_native_indices=list(range(7)), original_v3_component_proof_reused=True,
    component_interface=raw_component['component_interface'],
    raw_v4_objects=raw_v4_objects, raw_v4_direct_helper_joins=raw_direct_helper_joins,
    shared_value_math_object=shared_value_math_object, shared_value_interface=shared_value_interface,
    old_public_raw_v4_objects=old_public_raw_v4_objects,
    old_public_raw_v4_direct_helper_joins=old_public_direct_helper_joins,
    public_wrapper_byte_inverse=public_wrapper_byte_inverse,
    baseline_selector_function_joins=selector_joins,
    guard_shared_byte_authority=contract['guard_shared_byte_authority'],
    component_manifest=dict(path=contract['component_manifest'], **pin(manifest_path)),
    component_code_proof=component['members']['out/native-cpu-v4-arithmetic/code-proof.json'],
    legacy_math_component_manifest=dict(path=str(legacy['math_manifest'].relative_to(R)), **pin(legacy['math_manifest'])),
    legacy_optional_math_code_proof=legacy['math_code_proof'],
    checked_optional_math_splitter_retention_dependencies=36,
    expected_retained_input_functions=expected_functions,
    unchanged_callback_source=unchanged_sources, unchanged_bridge_headers=raw_component['unchanged_bridge_headers'],
    exact_two_header_inverse=True, admitted_guard_source_inverse=legacy['source_derivation']['guard_inverse_restores_AP_header'],
    actual_hidden_admitted_worker_calls=mesa_worker_calls,
    positive_job_preflight_branch=positive_branches, job_observation_layout=layout,
    callback_live_job_register=callback_job, fp_scope=fp_scope,
    eleven_input_merged_code=merge_association, eleven_input_final_code=part_association,
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
    scope='Static AU eleven-input composition only. All first7 AT objects, both Mesa headers, old public AS raw3 and GCC v4/v3/dispatcher roles remain exact. Only admitted v4 callbacks select the new scoped raw3 and shared9 value interface. All8241 exact PLT32/-4 call fields join actual unique hidden core definitions; private worker calls, genuine incoming job pointer/positive flag dominance and full FP save/set/restore are separately checked. Old pointer-bank/public and new typed-value interfaces are recorded separately. Dispatcher permits only the checked GOTPCRELX relaxation; other relocation masks keep the original bounds. Final native EVEX stays within guarded old/new v4 ranges. No provider loading, actual job execution, new numerical/older-hardware fallback execution, ABBA/frame/scientific/release acceptance or qualification transfer is established; genuine-job C8 remains the first numerical gate.')
with (P / 'linked-native-proof.json').open('x') as stream:
    stream.write(json.dumps(result, indent=2) + '\n')
print(json.dumps(dict(passed=True, native_objects=len(rows),
    linked_functions=sum(len(row['functions']) for row in rows), dispatcher_joins=len(dispatch_rows),
    complete_SPIRV_modules=len(modules), runtime_symbols_resolve=True)))
