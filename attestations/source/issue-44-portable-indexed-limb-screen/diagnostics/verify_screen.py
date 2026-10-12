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

import resource, subprocess, sys
sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
resource.setrlimit(resource.RLIMIT_FSIZE,(64*1024*1024,64*1024*1024))
stage=root/'bin/linux-gcc/portable-indexed-limb-accumulation'
tools=root/'bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools'
paths={'baseline':base/'retained_endpoint_portable_normal_sum.spv',
       'candidate':stage/'retained_endpoint_portable_normal_sum.spv'}
reports=[]
for producer,path in paths.items():
    subprocess.run([str(tools/'spirv-val'),'--target-env','vulkan1.2',str(path)],check=True)
    assembly=stage/(producer+'-verified.spvasm');assert not assembly.exists()
    subprocess.run([str(tools/'spirv-dis'),str(path),'-o',str(assembly)],check=True)
    text=assembly.read_text()
    bodies=[(m.group(1),m.group(0)) for m in re.finditer(r'(%\S+) = OpFunction .*?OpFunctionEnd',text,re.S)]
    accumulations=[(name,body) for name,body in bodies if re.fullmatch(r'%RTAccumulateTerm(?:_\d+)?',name)]
    assert len(accumulations)==2,(producer,'two original accumulator specializations')
    details=[]
    for name,body in accumulations:
        calls=re.findall(r'= OpFunctionCall %\S+ (%\S+)',body)
        sums=[v for v in calls if v.startswith('%RetainedSum32')]
        radii=[v for v in calls if v.startswith('%RPUpAdd')]
        bounds=[v for v in calls if v.startswith('%RPBounded')]
        assert len(sums)==1 and body.count('OpLoopMerge')==1,(producer,name)
        assert len(radii)==1 and len(bounds)==5,(producer,name)
        details.append({'function':name,'loop_merges':1,'sum_calls':sums,'radius_calls':radii,
                        'bound_calls':bounds,'complete_calls':calls,
                        'selection_merges':body.count('OpSelectionMerge'),
                        'selects':len(re.findall(r'\bOpSelect\b',body)),
                        'function_variables':[line.strip() for line in body.splitlines() if 'OpVariable' in line],
                        'access_chains':[line.strip() for line in body.splitlines() if 'OpAccessChain' in line],
                        'trip_tests':[line.strip() for line in body.splitlines() if 'OpULessThan' in line]})
    reports.append({'producer':producer,'identity_and_controls':module(path),
                    'specializations':details,'whole_module_loop_merges':text.count('OpLoopMerge'),
                    'whole_module_functions':len(bodies)})
a,b=reports
for key in ['capabilities','integer_types','barriers','execution_modes','descriptor_binding_layouts','no_contraction_operation_shapes']:
    assert a['identity_and_controls'][key]==b['identity_and_controls'][key],key
assert a['whole_module_loop_merges']==b['whole_module_loop_merges']
assert a['whole_module_functions']==b['whole_module_functions']
report={'passed':True,'reports':reports,'both_specializations_one_loop_one_sum_body':True,
        'single_radius_five_ordered_source_bounds':True,'actual_interface_modes_capabilities_integer_types_barriers_NoContraction_shapes_exact':True,
        'same_module_bytes':a['identity_and_controls']['sha256']==b['identity_and_controls']['sha256'],
        'scope':'One actual staged normal Endpoint product/interface screen; indexed lowering and private-frame/machine costs require subsequent inspection. No numerical execution, runtime benefit or qualification claim.'}
(work/'SPIRV-screen.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'passed':True,'baseline_bytes':a['identity_and_controls']['bytes'],
                  'candidate_bytes':b['identity_and_controls']['bytes'],
                  'specializations':len(b['specializations']),
                  'one_loop_one_sum_each':True,'all_structural_controls_exact':True}))
