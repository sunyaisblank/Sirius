"""Reuse the original called-code reader and add intergroup def-use checks."""
from pathlib import Path
import hashlib, json, re, struct
R = Path.cwd(); W = Path(__file__).resolve().parent

def intergroup(ir, v, event, folder):
    assert v['actual_entry']['entry_sysv_integer_arguments_2_3_4'] == [64, 1, 1]
    bindings = json.loads((folder / 'bindings.json').read_text())
    assert bindings['dispatches'][0]['groups_x'] == 1
    assert len(event['functions']) == 1 and event['functions'][0]['name'] == 'cs_variant'
    lines = ir.splitlines()
    definitions = dict(re.findall(r'^  (%[A-Za-z0-9_.]+) = (.*)$', ir, re.M))
    groups = [name for name, body in definitions.items() if re.search(r'add <8 x i32> .*<i32 0, i32 1, i32 2, i32 3, i32 4, i32 5, i32 6, i32 7>', body)]
    assert len(groups) == 1
    group, = groups
    def path(start, target):
        queue = [start]; parent = {start: None}
        for name in queue:
            if name == target:
                result = []
                while name is not None:
                    result.append(name + ' = ' + definitions.get(name, '<function argument>'))
                    name = parent[name]
                return result[::-1]
            for ref in re.findall(r'%[A-Za-z0-9_.]+', definitions.get(name, '')):
                if ref not in parent: parent[ref] = name; queue.append(ref)
        return []
    assert path(group, '%5'), 'Lane group IDs must derive from the original first-group-X argument'
    assert definitions['%owner_remaining'] == 'sub i32 %8, %5'
    assert definitions['%owner_tail'] == 'icmp ult i32 %owner_remaining, 8'
    assert definitions['%owner_partials'] == 'select i1 %owner_tail, i32 %owner_remaining, i32 0'
    assert 'splat (i32 -1)' in definitions['%mask.0']
    assert 'insertelement <8 x i32> %mask.1, i32 0, i32 %loop_counter2.0' in ir
    assert ir.index('%mask.0 = phi <8 x i32>') < ir.index('@llvm.masked.gather')
    calls = []
    for line in lines:
        if ' call ' not in line: continue
        m = re.search(r'@llvm.masked.(gather|scatter)\.[^(]+\((.*)\)', line)
        if not m: continue
        masks = re.findall(r'<8 x i1> (%[A-Za-z0-9_.]+)', m[2])
        assert len(masks) == 1, line
        mask_path = path(masks[0], '%mask.0')
        assert mask_path, ('Actual memory mask must retain pre-prologue tail-mask dependence', line)
        pointers = re.findall(r'<8 x ptr> (%[A-Za-z0-9_.]+)', m[2]); assert len(pointers) == 1
        address_path = path(pointers[0], group)
        calls.append({'operation': m[1], 'line': line.strip(), 'tail_mask_path': mask_path,
                      'group_address_path': address_path})
    gathers = [x for x in calls if x['operation'] == 'gather' and x['group_address_path']]
    scatters = [x for x in calls if x['operation'] == 'scatter' and x['group_address_path']]
    assert gathers and scatters
    return {'structural_premise_observed': True, 'original_block_arguments': [64, 1, 1],
        'guard_original_dispatch_groups_x': 1, 'lane_group_id_definition': group + ' = ' + definitions[group],
        'original_first_group_X_def_use': path(group, '%5'),
        'actual_row_address_and_mask_examples': [gathers[0], scatters[0]],
        'masked_memory_operations_with_tail_dependence': len(calls),
        'row_dependent_gathers': len(gathers), 'row_dependent_scatters': len(scatters),
        'review_limit': 'Mask def-use checks bind emitted lowering to the independently reviewed source mask composition. Guard exits before retained arithmetic; no active eight-row or numerical result inferred.'}

# Transport's original diagnostic used nm without an explicit nm prefreeze.
# Rejoin both selected symbols directly to the already frozen ELF, without
# repeating nm, GDB, compilation, dispatch or any workload.
folder = W / 'transport'; before = json.loads((folder / 'inputs-before.json').read_text())
exe = R / 'bin/linux-gcc/native-retained-ownership/short_input_jit_probe'
b = exe.read_bytes(); recorded = before[str(exe.relative_to(R))]
assert len(b) == recorded['bytes'] and hashlib.sha256(b).hexdigest() == recorded['sha256']
assert b[:6] == b'\x7fELF\x02\x01'
offset = struct.unpack_from('<Q', b, 40)[0]; size, count = struct.unpack_from('<HH', b, 58)
assert size == 64
sections = [struct.unpack_from('<IIQQQQIIQQ', b, offset + i * size) for i in range(count)]
names = []
for section in sections:
    if section[1] != 2: continue
    strings = sections[section[6]]; text = b[strings[4]:strings[4]+strings[5]]
    assert section[9] == 24
    for position in range(section[4], section[4]+section[5], 24):
        index, info, other, shndx, value, length = struct.unpack_from('<IBBHQQ', b, position)
        if shndx: names.append(text[index:].split(b'\0', 1)[0].decode())
selected = {}
for part in ['VulkanDevice8DispatchE', 'VulkanDevice10LoadKernelE']:
    matches = [name for name in names if part in name and '.cold' not in name]
    assert len(matches) == 1
    selected[part] = matches[0]
join = {'passed': True, 'frozen_executable': recorded, 'selected_symbols_from_ELF': selected,
    'nm_prefrozen_for_transport': False, 'capture_repeated': False,
    'scope': 'Independent exact symbol-name join against the original frozen ELF closes the diagnostic lookup dependency; it does not retroactively claim an nm tool prefreeze.'}
(folder / 'symbol-lookup-crosscheck.json').write_text(json.dumps(join, indent=2) + '\n')

basis = R / 'out/native-cpu-row-ownership/ordinary/read_code.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '1e5b5853a193cd0c3924a88dadca0e2960e7d5c24bd2c34274662e206131c7a2'
source = basis.read_text()
source = source.replace("(W.parent / 'compiler-result.json')", "(R / 'out/native-cpu-row-ownership/compiler-result.json')")
source = source.replace("(W.parent / 'emission-result.json')", "(R / 'out/native-cpu-row-ownership/emission-result.json')")
assert source.count('    rows[kind] = {') == 1
source = source.replace('    rows[kind] = {', '    extra = intergroup(ir, v, e, f)\n    rows[kind] = {**extra, ')
source = source.replace('Ordinary code remains8-wide SIMD with one active logical invocation; stack scratch is8 copies. Code is not smaller. No speed, memory-win, numerical, benefit, frame or release conclusion follows.',
    'Ordinary code has eight distinct X-group owners with masked row addresses and private lane storage. Guard activates only one short-input lane and exits before retained arithmetic. No speed, active eight-row numerical, benefit, frame or release conclusion follows.')
exec(compile(source, str(basis), 'exec'))
