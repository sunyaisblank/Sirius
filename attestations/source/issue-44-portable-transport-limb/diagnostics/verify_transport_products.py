"""Check complete producer identities and actual binary interface controls."""
from pathlib import Path
import hashlib
import json
import re
import struct

root = Path.cwd()
work = Path(__file__).resolve().parent
build = root / 'bin/linux-gcc/portable-transport-limb-accumulation'
base = root / 'bin/linux-gcc/src/sirius/backend/retained'
checked = root / 'attestations/software-vulkan/dcd1368/sealed-stage-seven-owner/files/bin/linux-gcc/sealed-stage-seven-owner/products'

def arrays(path):
    result = {}
    for size, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', path.read_text(), re.S):
        words = [int(v[:-1], 0) for v in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u', body)]
        assert len(words) == int(size)
        data = struct.pack('<%dI' % len(words), *words)
        result[name] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
    return result

def module(path):
    data = path.read_bytes()
    words = struct.unpack('<%dI' % (len(data) // 4), data)
    assert words[0] == 0x07230203
    instructions = []
    offset = 5
    while offset < len(words):
        size, opcode = words[offset] >> 16, words[offset] & 0xffff
        assert size and offset + size <= len(words)
        instructions.append((opcode, words[offset + 1:offset + size]))
        offset += size
    constants = {args[1]: args[2] for opcode, args in instructions if opcode == 43 and len(args) == 3}
    types = {args[0]: (opcode, args[1:]) for opcode, args in instructions if 19 <= opcode <= 39}
    decorations = {}
    members = {}
    for opcode, args in instructions:
        if opcode == 71:
            decorations.setdefault(args[0], []).append(tuple(args[1:]))
        elif opcode == 72:
            members.setdefault((args[0], args[1]), []).append(tuple(args[2:]))

    def signature(identifier):
        opcode, args = types[identifier]
        own = tuple(sorted(decorations.get(identifier, [])))
        if opcode in (21, 22):
            return (opcode, args, own)
        if opcode in (23, 24):
            return (opcode, signature(args[0]), args[1:], own)
        if opcode == 28:
            return (opcode, signature(args[0]), constants[args[1]], own)
        if opcode == 29:
            return (opcode, signature(args[0]), own)
        if opcode == 30:
            return (opcode, tuple((signature(member), tuple(sorted(members.get((identifier, index), []))))
                                 for index, member in enumerate(args)), own)
        if opcode == 32:
            return (opcode, args[0], signature(args[1]), own)
        raise AssertionError(('unexpected interface type', path, opcode))

    variables = {args[1]: (args[0], args[2]) for opcode, args in instructions if opcode == 59}
    layouts = []
    for identifier, values in decorations.items():
        binding = [v[1] for v in values if v[0] == 33]
        if not binding:
            continue
        descriptor = [v[1] for v in values if v[0] == 34]
        assert len(binding) == len(descriptor) == 1
        variable_type, storage = variables[identifier]
        qualifiers = tuple(sorted(v for v in values if v[0] not in (33, 34)))
        layouts.append((descriptor[0], binding[0], storage, qualifiers, signature(variable_type)))
    assert layouts
    definitions = {args[1]: (opcode, args[0]) for opcode, args in instructions
                   if len(args) >= 2 and not 19 <= opcode <= 39 and args[0] in types}
    contraction = []
    for identifier, values in decorations.items():
        for value in values:
            if value[0] == 42:
                opcode, result_type = definitions[identifier]
                contraction.append((opcode, signature(result_type)))
    barriers = sorted(tuple(constants[v] for v in args) for opcode, args in instructions if opcode == 224)
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
            'capabilities': [args for opcode, args in instructions if opcode == 17],
            'integer_types': sorted({args[1:] for opcode, args in instructions if opcode == 21}),
            'barriers': barriers, 'execution_modes': [args[1:] for opcode, args in instructions if opcode == 16],
            'descriptor_binding_layouts': sorted(layouts), 'no_contraction_operation_shapes': sorted(contraction)}

old = arrays(base / 'retained_kernels.h')
new = arrays(build / 'products/retained_kernels.h')
assert old.keys() == new.keys() and len(old) == 40
expected_arrays = {'kTransportPortableShader', 'kTransportPortableFp64Shader', 'kTransportPortableNormalSumShader'}
expected_modules = {'retained_transport_portable.spv', 'retained_transport_portable_fp64.spv',
                    'retained_transport_portable_normal_sum.spv', 'retained_transport_portable_stage7_owner.spv'}
changed_arrays = {name for name in old if old[name] != new[name]}
assert changed_arrays == expected_arrays, changed_arrays
assert len([name for name in old if 'Shader' not in name]) == 7
prior = root / 'attestations/software-vulkan/661e157/portable-limb-accumulation/files/bin/linux-gcc/portable-execution-attribution/limb-screen/products'
prior_arrays = arrays(prior / 'retained_kernels.h')
assert all(new[name] == prior_arrays[name] for name in expected_arrays)
records = []
for path in sorted((build / 'products').glob('*.spv')):
    original = base / path.name
    if not original.exists():
        original = checked / path.name
    assert original.is_file()
    a, b = module(original), module(path)
    for field in ('capabilities', 'integer_types', 'barriers', 'execution_modes', 'descriptor_binding_layouts', 'no_contraction_operation_shapes'):
        assert a[field] == b[field], (path, field)
    same = a['sha256'] == b['sha256']
    assert same == (path.name not in expected_modules), path
    measured_identity = None
    if path.name in expected_modules:
        measured_identity = module(prior / path.name)
        assert b == measured_identity, path
    records.append({'module': path.name, 'baseline_path': str(original.relative_to(root)),
                    'baseline': a, 'candidate': b, 'same_bytes': same,
                    'selected_Transport_matches_measured_661': measured_identity is not None,
                    'measured_661': measured_identity})
assert len(records) == 35
native = [record for record in records if '_portable' not in record['module']]
assert len(native) == 16 and all(record['same_bytes'] for record in native)
assert {record['module'] for record in records if not record['same_bytes']} == expected_modules
for name, value in (('transport-baseline-arrays.json', old), ('transport-candidate-arrays.json', new)):
    (work / name).write_text(json.dumps({'arrays': value}, indent=2) + '\n')
result = {'passed': True, 'complete_products': 35,
          'all7program_arrays_and16native_modules_exact': True,
          'other31products_baseline_exact': True,
          'selected4Transport_products_measured661_exact': True,
          'actual_capabilities_integer_types_barriers_modes_binding_layouts_NoContraction_shapes_exact': True,
          'arrays': {'same': [name for name in old if old[name] == new[name]],
                     'changed': sorted(changed_arrays)}, 'records': records,
          'scope': 'Complete producer products and structural controls with exact selected-module reuse. No new consumer execution, runtime benefit, full science or release claim.'}
(work / 'transport-generated-products-check.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'products': 35, 'unchanged_modules': 31, 'unchanged_arrays': 37,
                  'selected_Transport_matches_measured661': 4, 'all_structural_controls_exact': True}))
