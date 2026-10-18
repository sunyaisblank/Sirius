"""Passive saved-file/recorded-birth review; no project imports or external commands."""
from pathlib import Path
import hashlib, json, re

R = Path('/home/astra/.project/Sirius')
W = R / 'out/native-fp32-packed-soa-compiler-source'
P = R / 'out/native-fp32-packed-soa-publication'
A = R / 'attestations/qualification/5b2df58b/packed-retained-scratch-rejected/files/out/packed-retained-scratch-screen'
REV = '5b2df58bb9f24e7db5f7fbc0c38bdd832a58900e'
PATTERN = re.compile(r'^\s*(%ScratchAddress\S*) = OpFunction\b.*?OpFunctionEnd', re.M | re.S)


def identity(path):
    path = Path(path)
    data = path.read_bytes()
    return {'resolved_path': str(path.resolve(strict=True)), 'bytes': len(data),
            'sha256': hashlib.sha256(data).hexdigest()}


def pin(path):
    return {'path': str(path), **identity(path)}


def read(path):
    return json.loads(path.read_text())


def digest(text):
    return hashlib.sha256(text.encode()).hexdigest()


def address(path, registers, candidate):
    text = path.read_text()
    matches = list(PATTERN.finditer(text))
    assert len(matches) == 1
    found = matches[0]
    rows = [line.strip() for line in found.group(0).splitlines()]
    definitions = dict(re.findall(r'^\s*(%\S+) = (.+)$', text, re.M))
    tokens = [line.split() for line in rows]
    f, _, operation, u32, mode, signature = tokens[0]
    assert operation == 'OpFunction' and mode == 'None'
    assert definitions[u32] == 'OpTypeInt 32 0'
    assert definitions[signature].split() == ['OpTypeFunction', u32, u32, u32, u32, u32]
    parameters = []
    for row in tokens[1:5]:
        ident, equal, op, type_ = row
        assert equal == '=' and op == 'OpFunctionParameter' and type_ == u32
        parameters.append(ident)
    index, component, row_base, prefix = parameters
    def constant(ident, value):
        assert definitions[ident].split() == ['OpConstant', u32, str(value)]
    label = tokens[5][0]
    assert tokens[5] == [label, '=', 'OpLabel']
    temporary, _, _, pointer, storage = tokens[6]
    assert tokens[6][2] == 'OpVariable' and storage == 'Function'
    assert definitions[pointer].split() == ['OpTypePointer', 'Function', u32]
    if candidate:
        compare = tokens[7]
        condition, _, op, boolean, actual_component, four = compare
        assert op == 'OpIEqual' and actual_component == component
        assert definitions[boolean] == 'OpTypeBool'
        constant(four, 4)
        merge = tokens[8][1]
        yes, no = tokens[9][2:]
        assert tokens[8] == ['OpSelectionMerge', merge, 'None']
        assert tokens[9] == ['OpBranchConditional', condition, yes, no]
        assert tokens[10] == [yes, '=', 'OpLabel']
        assert tokens[11][:2] == ['OpStore', temporary]
        three = tokens[11][2]
        constant(three, 3)
        assert tokens[12] == ['OpBranch', merge]
        assert tokens[13] == [no, '=', 'OpLabel']
        assert tokens[14] == ['OpStore', temporary, component]
        assert tokens[15] == ['OpBranch', merge]
        assert tokens[16] == [merge, '=', 'OpLabel']
        plane = tokens[17][0]
        assert tokens[17] == [plane, '=', 'OpLoad', u32, temporary]
        multiplied, _, op, type_, source, stride = tokens[18]
        assert op == 'OpIMul' and type_ == u32 and source == plane
        constant(stride, registers)
        total = tokens[19][0]
        assert tokens[19] == [total, '=', 'OpIAdd', u32, multiplied, index]
        assert tokens[20] == ['OpReturnValue', total]
        assert tokens[21] == ['OpFunctionEnd'] and len(tokens) == 22
    else:
        # Explicit complete baseline graph, not assumed from old source.
        scaled, _, op, type_, source, four = tokens[7]
        assert op == 'OpIMul' and type_ == u32 and source == index
        constant(four, 4)
        condition, _, op, boolean, source, comparison_four = tokens[8]
        assert op == 'OpIEqual' and source == component and comparison_four == four
        assert definitions[boolean] == 'OpTypeBool'
        merge = tokens[9][1]
        yes, no = tokens[10][2:]
        assert tokens[9] == ['OpSelectionMerge', merge, 'None']
        assert tokens[10] == ['OpBranchConditional', condition, yes, no]
        assert tokens[11] == [yes, '=', 'OpLabel']
        assert tokens[12][:2] == ['OpStore', temporary]
        constant(tokens[12][2], 3)
        assert tokens[13] == ['OpBranch', merge]
        assert tokens[14] == [no, '=', 'OpLabel']
        assert tokens[15] == ['OpStore', temporary, component]
        assert tokens[16] == ['OpBranch', merge]
        assert tokens[17] == [merge, '=', 'OpLabel']
        plane = tokens[18][0]
        assert tokens[18] == [plane, '=', 'OpLoad', u32, temporary]
        total = tokens[19][0]
        assert tokens[19] == [total, '=', 'OpIAdd', u32, scaled, plane]
        assert tokens[20] == ['OpReturnValue', total]
        assert tokens[21] == ['OpFunctionEnd'] and len(tokens) == 22
    return text[:found.start()] + text[found.end():], found.group(0), {'function_id': f,
        'uint32_type': u32, 'function_type': signature, 'parameters': parameters,
        'function_local_address_temporary': temporary, 'blocks': [label, yes, no, merge],
        'comparison_result': condition, 'result': total, 'first_line': text.count('\n', 0, found.start()) + 1}


def main():
    target = P / 'independent-emitted-and-compiler-review.json'
    assert not target.exists()
    before = read(W / 'candidate/inputs-before.json')
    after = read(W / 'candidate/inputs-after.json')
    preflight = read(W / 'root-preflight.json')
    assert len(before) == len(after) == 4491 and before == after
    for key, row in before.items():
        actual = identity(Path(key) if Path(key).is_absolute() else R / key)
        assert actual == row, key
    assert len(preflight['files']) == 4490
    for key, row in preflight['files'].items():
        assert identity(key) == row, key
    owner = read(W / 'candidate/owner.json')
    assert owner['status'] == 'terminal' and owner['passed'] is True and owner['child_exit'] == 0
    assert owner['source'] == {'revision': REV, 'status': ''}
    assert owner['source_unchanged'] is True and owner['whole_inputs_unchanged'] is True
    assert owner['stop_reason'] is None and owner['cleanup_errors'] == [] and owner['remaining_owned_processes'] == []
    assert owner['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
    assert owner['declared_output_environment'] == {'GTEST_OUTPUT': ''} and owner['provider_mappings'] == []
    assert owner['boot_id'] == Path('/proc/sys/kernel/random/boot_id').read_text().strip()
    births = []
    for row in owner['observed_births']:
        assert row['absent'] is True
        path = Path('/proc') / str(row['pid']) / 'stat'
        try:
            stat = path.read_text()
            current = stat[stat.rfind(')') + 2:].split()[19]
        except FileNotFoundError:
            current = None
        assert str(current) != str(row['saved_start_ticks'])
        births.append({'pid': row['pid'], 'saved_start_ticks': str(row['saved_start_ticks']),
                       'fresh_current_start_ticks': current, 'original_birth_absent': True})
    assert len(births) == 6
    samples = [json.loads(line) for line in (W / 'candidate/samples.jsonl').read_text().splitlines()]
    assert max(row['rss_kib'] for row in samples) == owner['peak_sampled_rss_kib'] == 249992
    assert min(row['available_kib'] for row in samples) == owner['minimum_sampled_available_kib']
    recipe = read(W / 'recipes.json')
    commands = read(W / 'candidate/command-results.json')
    product = read(W / 'candidate/products.json')
    assert product['schema'] == 'sirius-native-packed-soa-five-compiler-products-v1'
    assert product['passed'] is True and product['runtime_admitted'] is False and product['emitted_address_semantics_admitted'] is False
    assert product['optimizer_invoked'] is False and product['actual_tool_calls'] == 25
    assert product['recipe'] == {'path': str(W / 'recipes.json'), **{k: identity(W / 'recipes.json')[k] for k in ['bytes','sha256']}}
    cases = []
    readable = []
    expected_commands = []
    tools = recipe['tools']
    source = Path(recipe['source_directory'])
    for row in recipe['profiles']:
        name = row['name']
        paths = {suffix: W / 'candidate' / (name + suffix) for suffix in recipe['expected_products_per_profile']}
        for suffix, path in paths.items():
            assert product['products'][name][suffix] == {'path': str(path), **{k: identity(path)[k] for k in ['bytes','sha256']}}
        for label, argv in [
          ('compile', [tools['compiler'],str(source/(name+'.slang')),*row['definitions'],'-I',str(source),'-O0','-target','spirv','-profile','spirv_1_5','-entry','ComputeMain','-stage','compute','-denorm-mode-fp32','preserve','-o',str(paths['.compiler.spv'])]),
          ('raw-disassemble',[tools['disassembler'],str(paths['.compiler.spv']),'-o',str(paths['.raw.spvasm'])]),
          ('assemble',[tools['assembler'],'--target-env','spv1.5',str(paths['.spvasm']),'-o',str(paths['.spv'])]),
          ('validate',[tools['validator'],'--target-env','vulkan1.2',str(paths['.spv'])]),
          ('final-disassemble',[tools['disassembler'],str(paths['.spv']),'-o',str(paths['.final.spvasm'])])]:
            expected_commands.append((name+'-'+label,argv))
        annotated = paths['.spvasm'].read_text()
        inverse = re.sub(r'^.*OpCapability RoundingModeRTE\n','',annotated,flags=re.M)
        inverse = re.sub(r'^.*OpExecutionMode %\S+ RoundingModeRTE 32\n','',inverse,flags=re.M)
        assert inverse == paths['.raw.spvasm'].read_text()
        old = A / 'posttools/disassemble' / ('candidate-'+name+'.spvasm')
        registers = row['general_table']['registers']
        left, old_address, old_ids = address(old, registers, False)
        right, new_address, new_ids = address(paths['.final.spvasm'], registers, True)
        assert left == right, ('outside-ScratchAddress', name)
        text = paths['.final.spvasm'].read_text()
        assert re.findall(r'OpExecutionMode %\S+ LocalSize (.+)', text) == ['64 1 1']
        declarations = dict(re.findall(r'^\s*(%\S+) = (.+)$', text, re.M))
        scratch = re.findall(r'OpName (%\S+) "retainedScratch"', text)
        assert len(scratch) == 1
        variable = declarations[scratch[0]].split()
        assert variable[0] == 'OpVariable' and variable[2] == 'Workgroup'
        pointer = declarations[variable[1]].split()
        array = declarations[pointer[2]].split()
        assert pointer[:2] == ['OpTypePointer','Workgroup'] and array[0] == 'OpTypeArray'
        assert declarations[array[1]] == 'OpTypeInt 32 0'
        length = declarations[array[2]].split()
        assert length[0] == 'OpConstant' and length[2] == str(4*registers)
        assert declarations[length[1]] == 'OpTypeInt 32 1'
        assert paths['.spv'].stat().st_size == (A/'candidate'/(name+'.spv')).stat().st_size
        cases.append({'name':name,'registers':registers,'scratch_uint_words':4*registers,
          'old_packed_spv':pin(A/'candidate'/(name+'.spv')),'candidate_spv':pin(paths['.spv']),
          'old_packed_asm':pin(old),'candidate_asm':pin(paths['.final.spvasm']),
          'outside_ScratchAddress_byte_exact_sha256':digest(left),'candidate_address':new_ids,
          'all_other_function_bodies_globals_types_constants_names_decorations_and_modes_byte_exact':True,
          'local_size':[64,1,1],'module_instruction_count':len([x for x in text.splitlines() if x.strip() and not x.lstrip().startswith(';')])})
        readable.append({'name':name,'old_address':old_address,'candidate_address':new_address,
          'old_ids':old_ids,'candidate_ids':new_ids,'outside_function_text_sha256':digest(left)})
    assert len(commands) == len(expected_commands) == 25
    for actual, (label, argv) in zip(commands,expected_commands):
        assert actual['label'] == label and actual['argv'] == argv and type(actual['exit_code']) is int and actual['exit_code'] == 0
        for stream in ['stdout','stderr']:
            path=Path(actual[stream]['path'])
            assert actual[stream] == {'path':str(path),**{k:identity(path)[k] for k in ['bytes','sha256']}}
    diff = P / 'source/address-graphs.json'
    diff.write_text(json.dumps(readable,indent=2)+'\n')
    result={'schema':'sirius-native-packed-soa-independent-emitted-and-compiler-review-v1',
      'passed':True,'source_revision':REV,'reviewer':'/root/sanitizer_gate_contract',
      'scope':'Independent actual terminal/owner/emitted-file review; reviewer authored compiler source, ROOT alone operated tools. No independent runtime/numerical/benefit/adoption admission.',
      'critical_pins':{key:pin(W/path) for key,path in [('preflight','root-preflight.json'),('actions','actions.json'),('recipes','recipes.json'),('compiler_owner','candidate/owner.json'),('before','candidate/inputs-before.json'),('after','candidate/inputs-after.json'),('command_results','candidate/command-results.json'),('products','candidate/products.json')]},
      'reader':pin(Path(__file__)),'address_graphs':pin(diff),
      'terminal':{'static_current_rows':4490,'before_after_current_rows_each':4491,'actual_commands':25,'wall_seconds':owner['wall_seconds'],'peak_sampled_rss_kib':owner['peak_sampled_rss_kib'],'minimum_sampled_available_kib':owner['minimum_sampled_available_kib'],'sample_count':len(samples),'fresh_recorded_births_absent':births,'cleanup_errors':[],'remaining_owned_processes':[]},
      'outside_helper_functions_alpha_exact':True,'outside_helper_text_byte_exact':True,'public_ABI_equal':True,'packed_component_major_address_validated':True,
      'comparison_strength':'For all5 pairs, entire saved assembly text after removing exactly one complete ScratchAddress function is byte-identical. Identity mapping of all remaining IDs/types/constants/function names is exact; no opaque AST equivalence or alpha-normalizer assumption. Includes ScratchRead/ScratchWrite/WriteValue/barriers/atomics/CFG, declarations/interfaces/modes. Complete old/new typed22-instruction address bodies checked, no unmatched instruction or side effect.',
      'address_contract':'uint32 plane=(component==4?3:component), product=plane*R, address=product+index; all4blocks/merge/control edges/local store/load exact checked. Under unchanged SOURCE preconditions component0..4/index<R, R459/612/815/362/629, address<4R<=3260: no uint32 overflow, bijective physical planes, logical3/4 same-register alias retained. rowBase/prefix parameters remain unused by this helper. No arbitrary malformed or stale payload equality promised.',
      'cases':cases,'decision':'EMITTED_AND_COMPILER_CONTRACT_PASS_RUNTIME_COST_UNPROVEN',
      'runtime_admitted':False,'production_adopted':False,
      'limits':['Source/component-major address validity composes with prior admitted logical validity0/1/error-sign/single-writer/barrier preconditions; not a whole-program formal/numerical proof.','All5 SPV sizes unchanged against packed AoS; no driver LDS/bank conflict/occupancy/cache/RSS/frame/runtime performance inference.','Only exact recorded Linux births rechecked; no global process-absence claim.','Other30 cached SPVs/40arrays/7tables/ownerheader composition remains future producer authority, not generated here. Original native controls+sole22gate ABBA still required before adoption.']}
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'receipt':pin(target),'address_graphs':pin(diff)},indent=2))


if __name__ == '__main__':
    main()
