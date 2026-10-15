"""Original full35/profile predicates; only six ASM reuse and signed typing adapt."""
from pathlib import Path
import hashlib, re, subprocess

def record(path):
    path=Path(path)
    with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
    return {'bytes':path.stat().st_size,'sha256':digest,'resolved_path':str(path.resolve(strict=True))}

def content(path):
    return {k:v for k,v in record(path).items() if k!='resolved_path'}

def estimated_divide_sites(stock,candidate,camera):
    from interface_signature import function_body
    if camera:
        import camera_profile
        return camera_profile.inspect(stock,candidate)
    import helper_normalize,graph
    helpers=helper_normalize.compare_helpers(stock,candidate)
    return {'inlined_camera':False,'complete_named_Divide_graph':graph.inspect(stock,candidate,function_body),
            'all_present_non_Divide_helpers_equal':helpers}


def semantic_instructions(text):
    ids=re.findall(r'^\s*(%\S+) =',text,re.M);assert len(ids)==len(set(ids))
    aliases={identifier:'id'+str(i) for i,identifier in enumerate(ids)}
    return [[aliases[token] if token.startswith('%') else token for token in re.findall(r'"(?:[^"\\]|\\.)*"|[^\s]+',line)]
            for line in text.splitlines() if line.strip() and not line.lstrip().startswith(';')]


def policy_receipt(module,raw,assembly,normal,camera):
    # The compiler owner executed the frozen policy. This reader joins its exact
    # persisted preimage/assembler/readback; no policy rerun.
    values=re.findall(r'^\s*(%\S+) = OpFDiv %\S+ .*$',raw,re.M)
    assert len(values)==(2 if camera else 1)
    annotation=''.join('\n               OpDecorate '+value+' NoContraction' for value in values)
    point=list(re.finditer(r'^.*OpDecorate .*?$',raw,re.M))[-1].end()
    expected=raw[:point]+annotation+raw[point:]
    if normal:
        entry=re.search(r'OpEntryPoint GLCompute (%\S+)',raw)[1]
        local=re.search(r'^.*OpExecutionMode '+re.escape(entry)+r' LocalSize.*$',raw,re.M)[0]
        expected=expected.replace('OpCapability Shader','OpCapability Shader\n               OpCapability RoundingModeRTE',1).replace(local,local+'\n               OpExecutionMode '+entry+' RoundingModeRTE 32',1)
    assert assembly==expected
    assert semantic_instructions(assembly)==semantic_instructions(module)
    actual_values=re.findall(r'^\s*(%\S+) = OpFDiv %\S+ .*$',module,re.M)
    assert len(actual_values)==len(values)
    assert all(re.search(r'OpDecorate '+re.escape(value)+r' NoContraction(?:\s|$)',module) for value in actual_values)
    return {'policy_completed_compiler_receipt_reused':True,'exact_annotation_inverse_count':len(values),
            'original_RTE_lines_added':2 if normal else 0,'complete_assembler_readback_instructions_equal':len(semantic_instructions(module))}


def products_review(P, R, W, T, original_work):
    from elf_arrays import arrays
    from interface_signature import interface_signature
    canonical = R / 'bin/linux-gcc/src/sirius/backend/retained'
    names = sorted(path.name for path in canonical.glob('*.spv'))
    assert len(names) == 35 and set(P['changed_modules']) <= set(names)
    rows = []
    for mode in ('stock', 'candidate'):
        assert sorted(path.name for path in (T / mode).glob('*.spv') if not path.name.endswith('.compiler.spv')) == names
        for path in (T / mode).iterdir():
            assert path.is_file() and not path.is_symlink()
    def barriers(text):
        constants = dict(re.findall(r'(%\S+) = OpConstant %\S+ (\d+)', text))
        return [[op, [constants[arg] for arg in args.split()]] for op, args in
                re.findall(r'^\s*(OpControlBarrier|OpMemoryBarrier) (.+)$', text, re.M)]
    for name in names:
        paths = {mode: T / mode / name for mode in ('stock', 'candidate')}
        assert paths['stock'].read_bytes() == (canonical / name).read_bytes()
        changed = paths['stock'].read_bytes() != paths['candidate'].read_bytes()
        assert changed == (name in P['changed_modules']), name
        texts = {}; assembly_records = {}
        for mode, path in paths.items():
            key = mode + '/' + name
            if key in P['saved_assembly']:
                target = Path(P['saved_assembly'][key]['path'])
                assert content(target) == P['saved_assembly'][key]['pin']
            else:
                target = W / 'analysis' / mode / (name + 'asm')
                target.parent.mkdir(parents=True, exist_ok=True)
                assert not target.exists()
                subprocess.run([P['disassembler'], str(path), '-o', str(target)], cwd=R, check=True)
            texts[mode] = target.read_text()
            assembly_records[mode] = {'assembly':record(target),'source_module':record(path),'reused_saved_assembly':key in P['saved_assembly']}
        stock, candidate = texts['stock'], texts['candidate']
        assert interface_signature(stock) == interface_signature(candidate), name
        modes = lambda text: sorted(re.findall(r'OpExecutionMode %\S+ (.*)$', text, re.M))
        assert modes(stock) == modes(candidate) and barriers(stock) == barriers(candidate), name
        footprint = None
        if changed:
            normal = name.endswith('_normal_sum.spv')
            caps = re.findall(r'OpCapability (\S+)', candidate)
            assert sorted(caps) == sorted(['Shader'] + (['RoundingModeRTE'] if normal else [])), name
            assert len(caps) == (2 if normal else 1)
            assert re.findall(r'OpTypeFloat (\d+)',candidate)==['32'],name
            assert set(re.findall(r'OpTypeInt (\d+) [01]',candidate)) == {'32'}, name
            assert not re.search(r'\b(?:RelaxedPrecision|FPFastMathMode|FPFastMathDefault|OpFmaKHR|OpConvertFTo[SU]|OpConvert[SU]ToF|NoUnsignedWrap|NoSignedWrap)\b',candidate)
            camera=name.startswith(('retained_camera_','retained_ray_camera_'))
            persisted_policy=policy_receipt(candidate,(T/'candidate'/name).with_suffix('.raw.spvasm').read_text(),(T/'candidate'/name).with_suffix('.spvasm').read_text(),normal,camera)
            footprint=estimated_divide_sites(stock,candidate,camera)
            footprint['persisted_compiler_policy_join']=persisted_policy
        rows.append({'module': name, 'changed': changed, 'stock': record(paths['stock']), 'candidate': record(paths['candidate']),
                     'complete_public_shared_type_layout_modes_and_ordered_barriers_equal': True, 'divide_footprint': footprint, 'assembly_binary_joins':assembly_records})
    headers = {}
    for name, count in [('retained_kernels.h', 40), ('retained_owner_kernels.h', 2)]:
        expected = arrays(canonical / name)
        stock = arrays(T / 'stock' / name); candidate = arrays(T / 'candidate' / name)
        assert len(stock) == count and stock == expected and stock.keys() == candidate.keys()
        changed = [key for key in stock if stock[key] != candidate[key]]
        assert set(changed) == (set(P['changed_arrays']) if count == 40 else set())
        assert all(stock[key] == candidate[key] for key in stock if key.endswith('Program'))
        strip = lambda text: re.sub(r'inline constexpr std::array<std::uint32_t, \d+> (\w+)\{\{.*?\}\};', r'ARRAY \1', text, flags=re.S)
        assert strip((canonical / name).read_text()) == strip((T / 'stock' / name).read_text()) == strip((T / 'candidate' / name).read_text())
        headers[name] = {'changed_arrays': changed, 'stock': record(T / 'stock' / name), 'candidate': record(T / 'candidate' / name)}
    for mode in ('stock','candidate'):
        inventory={'main':arrays(T/mode/'retained_kernels.h'),'owner':arrays(T/mode/'retained_owner_kernels.h')}
        for module,row in P['module_arrays'].items():
            assert inventory[row['header']][row['array']]==(T/mode/module).read_bytes(),(mode,module)
    for name in ('retained_owner_layout.h', 'retained_owner_layout.slang'):
        assert (canonical / name).read_bytes() == (T / 'stock' / name).read_bytes() == (T / 'candidate' / name).read_bytes()
    return {'passed': True, 'stock_reuse':record(W/'stock-reuse-acceptance.json'), 'modules': rows, 'headers': headers, 'changed_modules': 17, 'unchanged_modules': 18,
            'main_arrays_changed': 17, 'main_arrays_unchanged': 23, 'owner_arrays_unchanged': 2,
            'all_seven_tables_and_nonarray_ABI_constants_equal': True,
            'scope': 'Complete original stock35+new candidate35/40+2/7tables source/profile evidence; separately owned compiler remains truthful if this reader refuses. No integrated numerical or benefit verdict.'}
