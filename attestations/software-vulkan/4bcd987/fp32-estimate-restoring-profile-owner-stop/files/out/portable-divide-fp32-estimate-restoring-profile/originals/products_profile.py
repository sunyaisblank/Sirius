"""Original full35/profile predicates; only six ASM reuse and signed typing adapt."""
from pathlib import Path
import hashlib, re, subprocess

def record(path):
    path=Path(path)
    with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
    return {'bytes':path.stat().st_size,'sha256':digest,'resolved_path':str(path.resolve(strict=True))}

def content(path):
    return {k:v for k,v in record(path).items() if k!='resolved_path'}

def estimated_divide_sites(stock, candidate, camera):
    from interface_signature import function_body, normalized_function
    if camera:
        import camera_profile
        return camera_profile.inspect(stock, candidate)
    import estimate_graph
    for name in ('PB32Multiply','PB32ProductResidual','PB32Unpack','PB32Round'):
        assert normalized_function(stock, name) == normalized_function(candidate, name), name
    return {'inlined_camera':False,'complete_named_Divide_graph':estimate_graph.inspect(stock,candidate,function_body),
            'original_multiply_residual_unpack_round_normalized_bodies_equal':True}

def products_review(P, R, W, T, original_work):
    from elf_arrays import arrays
    from interface_signature import interface_signature
    canonical = R / 'bin/linux-gcc/src/sirius/backend/retained'
    names = sorted(path.name for path in canonical.glob('*.spv'))
    assert len(names) == 35 and set(P['changed_modules']) <= set(names)
    rows = []
    for mode in ('stock', 'candidate'):
        assert sorted(path.name for path in (T / mode).glob('*.spv')) == names
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
            assert sorted(caps) == sorted(['Shader', 'Float64'] + (['RoundingModeRTE'] if normal else [])), name
            assert len(caps) == (3 if normal else 2)
            assert sorted(re.findall(r'OpTypeFloat (\d+)', candidate)) == (['32','64'] if normal else ['64']), name
            assert set(re.findall(r'OpTypeInt (\d+) [01]',candidate)) == {'32'}, name
            assert not re.search(r'\b(?:RelaxedPrecision|FPFastMathMode|FPFastMathDefault|OpFmaKHR|OpConvertFTo[SU])\b', candidate)
            footprint = estimated_divide_sites(stock, candidate, name.startswith(('retained_camera_', 'retained_ray_camera_')))
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
    for name in ('retained_owner_layout.h', 'retained_owner_layout.slang'):
        assert (canonical / name).read_bytes() == (T / 'stock' / name).read_bytes() == (T / 'candidate' / name).read_bytes()
    return {'passed': True, 'stock_reuse': record(original_work / 'stock-reuse-acceptance.json'), 'modules': rows, 'headers': headers, 'changed_modules': 17, 'unchanged_modules': 18,
            'main_arrays_changed': 17, 'main_arrays_unchanged': 23, 'owner_arrays_unchanged': 2,
            'all_seven_tables_and_nonarray_ABI_constants_equal': True,
            'scope': 'Complete preserved compiler products and fresh parse-only35module/profile/ABI evidence; old wrapper remains failed. No integrated numeric or benefit verdict.'}
