"""Bounded original default-profile products, producer and eight-control owner actions."""
from pathlib import Path
import ast, hashlib, json, re, shlex, shutil, subprocess, sys
sys.dont_write_bytecode = True
W = Path(__file__).resolve().parent
R = W.parents[1]
T = R / 'bin/linux-gcc/portable-divide-fp64-estimate-controls'
P = json.loads((W / 'packet.json').read_text())

def canonical_json(value):
    normalized = json.loads(json.dumps(value, allow_nan=False))
    return json.dumps(normalized, sort_keys=True, separators=(',', ':'), allow_nan=False)

def record(path):
    path = Path(path)
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': path.stat().st_size, 'sha256': digest,
            'resolved_path': str(path.resolve(strict=True))}

def content(path):
    return {k: v for k, v in record(path).items() if k != 'resolved_path'}

def write(path, value):
    path = Path(path)
    assert not path.exists(), path
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')

def verify_static():
    for path, pin in json.loads((W / 'source-freeze.json').read_text())['pins'].items():
        assert record(path) == pin, path

def accepted_folder(folder):
    folder = Path(folder)
    owner = json.loads((folder / 'owner.json').read_text())
    assert owner['status'] == 'terminal' and owner['passed'] is True and owner['child_exit'] == 0
    assert owner['stop_reason'] is None and owner['cleanup_errors'] == []
    assert owner['remaining_owned_processes'] == []
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
    return owner

def accepted(label):
    return accepted_folder(W / label)

def relocated(path):
    return Path(P['historical_host_archive']) / 'files' / Path(path).relative_to(R)

def verify_host_reuse():
    history = Path(P['historical_host_root'])
    accepted_folder(history / 'compile-device')
    seal = json.loads((history / 'compile-device/output-seal.json').read_text())
    for key, relative in [('source_freeze', 'source-freeze.json'),
                          ('operator_context', 'operator-context.json'),
                          ('inputs_before', 'compile-device/inputs-before.json')]:
        assert content(history / relative) == {k: seal[key][k] for k in ('bytes', 'sha256')}
    for original, pin in seal['outputs'].items():
        assert content(relocated(original)) == {k: pin[k] for k in ('bytes', 'sha256')}
    assert str(P['device_object']) in [str(relocated(name)) for name in seal['outputs']]
    pair = json.loads(Path(P['historical_pair_review']).read_text())
    assert pair['passed'] and pair['complete_55_enabled_fields_queue_extensions_order_and_feature_chain_equal']
    features = []
    for name in ('stock', 'candidate'):
        folder = history / ('capture-' + name + '-endpoint-normal')
        accepted_folder(folder)
        features.append(json.loads((folder / 'feature-review.json').read_text()))
        assert features[-1]['actual_public_create_arguments_enable_shaderInt64']
        assert features[-1]['successful_logical_device_used_for_actual_complete_shader']
    assert features[0] == features[1] == pair['normalized_feature_report']
    assert len(features[0]['features']) == 55 and features[0]['features']['shaderInt64'] == 1
    assert features[0]['features']['shaderFloat64'] == 1
    return {'device_object': record(P['device_object']), 'original_output_seal': record(history / 'compile-device/output-seal.json'),
            'historical_pair_review': record(P['historical_pair_review']),
            'scope': 'Exact accepted compiled Int64 device object/source and historical public-argument witness reused. No fresh numeric-process argument capture is claimed.'}

def verify_output_seal(label):
    accepted(label)
    seal = json.loads((W / label / 'output-seal.json').read_text())
    assert seal['label'] == label and seal['recipe'] == P['recipes'][label]
    assert seal['source_freeze'] == record(W / 'source-freeze.json')
    assert seal['operator_context'] == record(W / 'operator-context.json')
    assert seal['inputs_before'] == record(W / label / 'inputs-before.json')
    assert seal['outputs'] == {name: record(name) for name in P['outputs'][label]}
    return record(W / label / 'output-seal.json')

def join_runtime_inputs(label, owner):
    action = json.loads((W / 'actions.json').read_text())[label]
    folder = W / label
    assert owner['argv'] == [arg.replace('{folder}', str(folder)) for arg in action['argv']]
    before = json.loads((folder / 'inputs-before.json').read_text())
    for name in action['frozen_inputs']:
        path = Path(name)
        key = str(path.relative_to(R)) if path.is_relative_to(R) else str(path)
        assert before[key] == record(path), ('numerical consumed input differs from current sealed input', label, name)

def kernel_source_profile():
    old = Path(P['original_stock_work']) / 'stock/src/sirius/kernels'
    private = W / 'stock/src/sirius/kernels'
    canonical = R / 'src/sirius/kernels'
    names = sorted(p.name for p in old.iterdir())
    assert len(names) == 33 and names == sorted(p.name for p in private.iterdir())
    assert names == sorted(p.name for p in (W / 'candidate/src/sirius/kernels').iterdir())
    rows, deltas = {}, []
    for name in names:
        stock = private / name; current = canonical / name
        assert stock.read_bytes() == (old / name).read_bytes(), name
        if name != 'portable_binary32.h':
            assert stock.read_bytes() == (W / 'candidate/src/sirius/kernels' / name).read_bytes(), name
        different = stock.read_bytes() != current.read_bytes()
        if different:
            deltas.append(name)
            expected = P['known_outside_retained_source_deltas'][name]
            assert content(stock) == expected['original_private'] and content(current) == expected['canonical4498'], name
        rows[name] = {'original':record(old/name), 'private':record(stock), 'canonical4498':record(current), 'current_equal':not different}
    assert sorted(deltas) == sorted(P['known_outside_retained_source_deltas']) == ['gr_ellis_topology.slang','gr_integrator.slang','parity_probe.slang','trace.slang']
    # Literal Slang import/include closure only. Every directive is recognized;
    # macro/dynamic includes, include_next and alternate search names refuse.
    seen, pending, generated = set(), list(P['retained_source_roots']), set()
    while pending:
        name = pending.pop()
        if name in seen: continue
        seen.add(name); assert name in names
        text = (private/name).read_text()
        if name.endswith('.py'):
            tree = ast.parse(text)
            imports = []
            for node in ast.walk(tree):
                if isinstance(node,ast.Import): imports.extend(x.name for x in node.names)
                elif isinstance(node,ast.ImportFrom):
                    assert node.level == 0; imports.append(node.module)
                elif isinstance(node,ast.Call):
                    assert not (isinstance(node.func,ast.Name) and node.func.id in ('open','eval','exec','compile','__import__','getattr'))
                    assert not (isinstance(node.func,ast.Attribute) and node.func.attr in ('read_text','read_bytes','open','read','write_text','write_bytes','exec_module','import_module'))
            assert sorted(imports) == ['fractions','struct']
            continue
        text = re.sub(r'/\*.*?\*/|//[^\n]*','',text,flags=re.S)
        for line in text.splitlines():
            if re.search(r'\bimport\b',line):
                match = re.fullmatch(r'\s*import (\w+);\s*',line); assert match, line
                dependency = match[1]+'.slang'
                if dependency == 'retained_owner_layout.slang':
                    assert name == 'retained_owner.slang'; generated.add(dependency)
                else: pending.append(dependency)
            if re.search(r'^\s*#\s*(?:include|include_next|import)\b',line):
                if name == 'portable_binary32.h' and line in ('#include <bit>','#include <cstdint>'):
                    cpp_branch = text.split('#if defined(__cplusplus)',1)[1].split('#else',1)[0]
                    assert line in cpp_branch.splitlines()
                    continue  # Fixed C++-only headers, inactive for retained Slang.
                match = re.fullmatch(r'\s*#include "([A-Za-z0-9_.]+)"\s*',line); assert match, line
                pending.append(match[1])
    assert sorted(seen) == P['retained_source_closure'] and len(seen) == 14
    assert generated == {'retained_owner_layout.slang'} and not set(deltas)&seen
    assert all(rows[name]['current_equal'] for name in seen)
    builder = (W / 'stock/scripts/build-retained-kernels.py').read_text()
    assert builder == (R / 'scripts/build-retained-kernels.py').read_text()
    assert builder == (Path(P['original_stock_work']) / 'stock/scripts/build-retained-kernels.py').read_text()
    assert builder.count('spec_from_file_location("retained_program", source / "retained_program.py")') == 1
    assert 'owner_include = ["-I", str(destination.parent)] if stage_seven_owner else []' in builder
    assert '[compiler, str(source), *definitions, "-I", str(source.parent), *owner_include, optimization,' in builder
    # Both search roots are fixed by this source: private source directory and
    # owned generated owner-layout directory. No CLI source/include argument.
    cmake = (R / 'src/sirius/backend/CMakeLists.txt').read_text()
    deps = cmake.split('DEPENDS "${CMAKE_SOURCE_DIR}/scripts/build-retained-kernels.py"',1)[1].split('COMMENT "Building bounded retained camera and joint transport stages"',1)[0]
    assert sorted(re.findall(r'"\$\{SIRIUS_KERNEL_DIR\}/([^"/]+)"',deps)) == sorted(seen)
    layout = lambda source: source.split('def generate_owner_layout(',1)[1].split('\ndef compile_shader(',1)[0]
    assert layout(builder) == layout((W / 'candidate/scripts/build-retained-kernels.py').read_text())
    return {'historical33_original_private_sources':rows, 'complete_retained_import_include_closure':sorted(seen),
            'generated_owner_layout':sorted(generated),'inactive_cpp_only_standard_headers':['bit','cstdint'],'all14_consumed_source_files_equal_canonical4498':True,
            'current_equal_source_count':len(names)-len(deltas), 'exact4_current_deltas_outside_retained_closure':deltas,
            'builder':record(W/'stock/scripts/build-retained-kernels.py'),'CMake_dependency_authority':record(R/'src/sirius/backend/CMakeLists.txt'),
            'scope':'Historical4bcd complete stock sourcebase preserved, complete retained actual-read14 closure equal current4498. Four unrelated floor-repaired trace/probe files are frozen current deltas outside that closure. No current full trace/science qualification.'}

def verify_arithmetic():
    archive = Path(P['arithmetic_archive'])
    assert re.fullmatch('[0-9a-f]{40}', P['arithmetic_commit']), 'immutable arithmetic prerequisite pending'
    assert re.fullmatch('[0-9a-f]{64}', P['arithmetic_manifest_sha256']), 'immutable arithmetic manifest pending'
    assert content(archive / 'manifest.json')['sha256'] == P['arithmetic_manifest_sha256']
    manifest = json.loads((archive / 'manifest.json').read_text())
    for name, pin in manifest['files'].items():
        assert not (archive / name).is_symlink() and content(archive / name) == pin, name
    work = archive / 'files/out/portable-divide-fp64-estimate-arithmetic'
    historical = R / 'out/portable-divide-fp64-estimate-arithmetic'
    def join(name, pin):
        original = Path(name)
        if not original.is_absolute(): original = R / original
        saved = archive / 'files' / original.relative_to(R) if original.is_relative_to(R) else None
        actual = saved if saved is not None and saved.is_file() else original
        assert content(actual) == {key:pin[key] for key in ('bytes','sha256')}, name
        if actual == original: assert record(actual) == pin, name
        else: assert not actual.is_symlink() and pin['resolved_path'] == str(original), name
        return actual
    frozen = json.loads((work / 'source-freeze.json').read_text())['pins']
    assert len(frozen) == 1886
    for name,pin in frozen.items(): join(name,pin)
    result_path = work / 'arithmetic-result.json'
    assert content(result_path)['sha256'] == P['arithmetic_inner_result_sha256']
    result = json.loads(result_path.read_text())
    assert result['schema'] == 'sirius-fp64-estimate-original-arithmetic-result-v1'
    assert result['passed'] and result['original_records'] == 174268 and result['original_operation_count'] == 14
    assert result['complete_output_words_per_provider'] == 697072
    assert result['fresh_candidate_host_and_full_gpu_plus_reused_stock_equal_original_expected']
    assert result['stock_result_reused'] and result['candidate_macro'] == P['candidate_macro'] == 'SIRIUS_PB32_DIVIDE_FP64_ESTIMATE=1'
    assert result['original_corpus']['division_branch_coverage']['optional_finite_nonzero_branch_records'] == 12263
    assert (work / 'candidate/src/sirius/kernels/portable_binary32.h').read_bytes() == (W / 'candidate/src/sirius/kernels/portable_binary32.h').read_bytes()
    terminal = work / 'independent-terminal-review.json'
    assert content(terminal)['sha256'] == P['arithmetic_terminal_review_sha256'] and json.loads(terminal.read_text())['passed']
    assert len(result['owners']) == 10
    for label,pin in result['owners'].items():
        accepted_folder(work / label)
        assert content(work / label / 'owner.json') == {key:pin[key] for key in ('bytes','sha256')}
        before = json.loads((work / label / 'inputs-before.json').read_text())
        for name,pin in before.items(): join(name,pin)
    assert (work / 'execute-gpu-candidate/actual.bin').read_bytes() == (work / 'expected.bin').read_bytes()
    assert len((work / 'expected.bin').read_bytes()) == 697072 * 4
    source = (W / 'candidate/src/sirius/kernels/portable_binary32.h').read_bytes()
    insertion = (W / 'candidate-insertion.txt').read_bytes()
    assert source.count(insertion) == 1 and source.replace(insertion,b'',1) == (W / 'stock/src/sirius/kernels/portable_binary32.h').read_bytes()
    source_profile = kernel_source_profile()
    return {'result':record(result_path),'terminal_review':record(terminal),'manifest':record(archive/'manifest.json'),
            'all1886_original_static_and10terminal_owners_input_joins':True, 'retained_source_profile':source_profile,
            'scope':'Accepted original174268/14-op primitive corpus, exact historical stock reuse and fresh candidate host/provider output. No estimate hit-rate, retained-mode, integrated or benefit acceptance inferred.'}

def pin_context():
    verify_static(); verify_host_reuse(); verify_arithmetic()
    assert not T.exists() and not (W / 'operator-context.json').exists()
    assert not subprocess.check_output(['git', 'status', '--porcelain=v1'], cwd=R)
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=R, text=True).strip()
    assert (W / 'actions.json').read_bytes() == (W / 'actions.unpinned.json').read_bytes()
    assert revision == P['runner_preparation_revision']
    actions = json.loads((W / 'actions.json').read_text())
    assert all(action['revision'] == revision for action in actions.values())
    write(W / 'operator-context.json', {'revision': revision, 'actions': record(W / 'actions.json'),
          'source_freeze': record(W / 'source-freeze.json'), 'preparation': record(W / 'preparation.json')})
    print('Pinned original Float64-estimate profile; no compiler/provider launched.', flush=True)

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

def reuse_stock():
    archive = Path(P['original_stock_archive'])
    manifest_path = archive / 'manifest.json'
    assert content(manifest_path)['sha256'] == P['original_stock_manifest_sha256']
    manifest = json.loads(manifest_path.read_text())
    old_work = Path(P['original_stock_work'])
    accepted_folder(old_work / 'compile-products-stock')
    before = json.loads((old_work / 'compile-products-stock/inputs-before.json').read_text())
    tools = {}
    for name in P['stock_reuse_tools']:
        path = Path(name)
        old = next(value for value in before.values() if value['resolved_path'] == str(path.resolve(strict=True)))
        assert content(path) == {key: old[key] for key in ('bytes', 'sha256')}
        tools[name] = record(path)
    kernel_sources = {}
    original = old_work / 'stock/src/sirius/kernels'
    names = sorted(path.name for path in original.iterdir())
    assert len(names) == 33 and names == sorted(path.name for path in (W / 'stock/src/sirius/kernels').iterdir())
    for name in names:
        assert (original/name).read_bytes() == (W/'stock/src/sirius/kernels'/name).read_bytes()
        kernel_sources[name] = record(original / name)
    builders = [old_work / 'stock/scripts/build-retained-kernels.py', W / 'stock/scripts/build-retained-kernels.py', R / 'scripts/build-retained-kernels.py']
    assert len({content(path)['sha256'] for path in builders}) == 1
    source = Path(P['stock_reuse_directory'])
    target = T / 'stock'; assert not target.exists()
    products = [Path(name).name for name in P['outputs']['reuse-products-stock'] if Path(name).parent == target]
    assert len(products) == 39 and set(products) == {path.name for path in source.iterdir()}
    target.mkdir(parents=True)
    joined = {}
    for name in products:
        origin = source / name
        key = str(origin.relative_to(archive))
        assert content(origin) == manifest['files'][key]
        assert origin.read_bytes() == (R / 'bin/linux-gcc/src/sirius/backend/retained' / name).read_bytes()
        shutil.copyfile(origin, target / name)
        assert (target / name).read_bytes() == origin.read_bytes()
        joined[name] = {'source': record(origin), 'copy': record(target / name)}
    return {'passed': True, 'reused_files': 39, 'original_stock_owner': record(old_work / 'compile-products-stock/owner.json'),
            'original_stock_before': record(old_work / 'compile-products-stock/inputs-before.json'),
            'original_stock_after': record(old_work / 'compile-products-stock/inputs-after.json'),
            'original_manifest': record(manifest_path), 'kernel_sources': kernel_sources,
            'builder': record(builders[0]), 'five_tools': tools, 'products': joined, 'retained_source_profile':kernel_source_profile(),
            'preparation_cost_evidence_seconds': {'historical_stock35': 35.46063927599971, 'historical_candidate35': 34.996513119000156},
            'scope': 'Exact accepted original stock35modules+4headers/maps reused; no fresh stock compiler or runtime result. Full candidate35 builder retained to avoid a selective compiler/cache adapter; no future timing promise or benefit claim.'}

def products_review():
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
        texts = {}
        for mode, path in paths.items():
            target = T / 'analysis' / mode / (name + 'asm')
            target.parent.mkdir(parents=True, exist_ok=True)
            assert not target.exists()
            subprocess.run([P['disassembler'], str(path), '-o', str(target)], cwd=R, check=True)
            texts[mode] = target.read_text()
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
                     'complete_public_shared_type_layout_modes_and_ordered_barriers_equal': True, 'divide_footprint': footprint})
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
    return {'passed': True, 'stock_reuse': record(W / 'stock-reuse-acceptance.json'), 'modules': rows, 'headers': headers, 'changed_modules': 17, 'unchanged_modules': 18,
            'main_arrays_changed': 17, 'main_arrays_unchanged': 23, 'owner_arrays_unchanged': 2,
            'all_seven_tables_and_nonarray_ABI_constants_equal': True,
            'scope': 'Complete compiler/profile/ABI evidence only; no integrated numeric or benefit verdict.'}

def prepare_producer():
    assert json.loads((W / 'products-acceptance.json').read_text())['passed']
    verify_arithmetic(); verify_host_reuse()
    directory = T / 'controls'; assert not directory.exists()
    (directory / 'headers').mkdir(parents=True)
    for name in ('retained_kernels.h', 'retained_owner_kernels.h', 'retained_owner_layout.h'):
        shutil.copyfile(T / 'candidate' / name, directory / 'headers' / name)
    shutil.copyfile(P['baseline_origin'], directory / 'original_baseline_controls')
    (directory / 'original_baseline_controls').chmod(0o755)
    for relative, source in P['resource_sources'].items():
        target = directory / 'resources' / relative
        target.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(source, target)
    protected = json.loads((W / 'protected-before.json').read_text())
    for name, pin in protected.items(): assert content(R / name) == pin, name

def producer_feature_review():
    needed = json.loads(Path(P['original_symbols']).read_text())['required_strong_sirius_definitions']
    strong = lambda text: {row.split()[2] for row in text.splitlines() if len(row.split()) == 3 and row.split()[1] in ('T', 'D', 'B', 'R')}
    symbols = {}
    for name, path in [('device-object', Path(P['device_object'])), ('candidate', T / 'controls/original_candidate_controls')]:
        text = subprocess.check_output(['/usr/bin/nm', '-g', '--defined-only', str(path)], text=True)
        assert set(needed) <= strong(text); symbols[name] = text
    linkmap = (T / 'controls/candidate.map').read_text()
    assert 'libsirius_backend.a(vulkan_device.cpp.o)' not in linkmap
    assert 'libsirius_backend.a(retained_compute.cpp.o)' not in linkmap
    assert str(R / 'bin/linux-gcc/src/sirius/backend/CMakeFiles/sirius_backend.dir/vulkan/vulkan_device.cpp.o') not in linkmap
    argv = P['recipes']['link-candidate'][0]
    assert argv.count(P['device_object']) == 1
    assert argv.index(P['device_object']) < argv.index(str(R / 'bin/linux-gcc/src/sirius/backend/libsirius_backend.a'))
    assert not any('--whole-archive' in arg for arg in argv)
    producer = json.loads((W / 'producer-acceptance.json').read_text())
    assert producer['passed'] and set(producer['changed_arrays']) == set(P['changed_arrays'])
    assert producer['unchanged_array_count'] == 23
    for name in ('kTransportShader', 'kEndpointShader'):
        assert producer['actual_ELF_owner_arrays']['baseline'][name]['sha256'] == producer['actual_ELF_owner_arrays']['candidate'][name]['sha256']
    return {'passed': True, 'all28device_definitions_in_object_and_actual_candidate': True,
            'old_vulkan_and_compute_archive_members_absent': True, 'link_uses_default_archive_extraction': True,
            'feature_host': verify_host_reuse(), 'nm_outputs': symbols,
            'scope': 'Candidate linked to identical admitted diagnostic device object. Prior public Float64 arguments remain historical; the reused diagnostic object also enables unused Int64. Candidate shaders require noInt64. Original baseline ELF remains unchanged; public optional-provider fallback/admission is not established.'}

def ancestor_labels(label):
    ordered = []
    seen = set()
    def visit(current):
        for parent in P['dependencies'][current]:
            if parent not in seen:
                seen.add(parent); visit(parent); ordered.append(parent)
    visit(label)
    return ordered

def phase(label):
    verify_static()
    context = json.loads((W / 'operator-context.json').read_text())
    assert context['actions'] == record(W / 'actions.json') and context['preparation'] == record(W / 'preparation.json')
    assert context['source_freeze'] == record(W / 'source-freeze.json')
    for prior in ancestor_labels(label): verify_output_seal(prior)
    assert all(not Path(path).exists() for path in P['outputs'][label])
    if label == 'reuse-products-stock': write(W / 'stock-reuse-acceptance.json', reuse_stock())
    elif label == 'verify-products': write(W / 'products-acceptance.json', products_review())
    elif label == 'prepare-producer': prepare_producer()
    else:
        for argv in P['recipes'][label]: subprocess.run(argv, cwd=R, check=True)
        if label == 'verify-producer': write(W / 'producer-feature-acceptance.json', producer_feature_review())
    write(W / label / 'output-seal.json', {'label': label, 'recipe': P['recipes'][label],
          'source_freeze': record(W / 'source-freeze.json'), 'operator_context': record(W / 'operator-context.json'),
          'inputs_before': record(W / label / 'inputs-before.json'),
          'outputs': {name: record(name) for name in P['outputs'][label]}})

def pin_runtime():
    verify_static(); verify_arithmetic(); verify_host_reuse()
    assert not (W / 'runtime-context.json').exists()
    for label in P['producer_labels']: verify_output_seal(label)
    assert json.loads((W / 'producer-acceptance.json').read_text())['passed']
    assert json.loads((W / 'producer-feature-acceptance.json').read_text())['passed']
    write(W / 'runtime-context.json', {'profile': record(W / 'profile.json'),
          'producer': record(W / 'producer-acceptance.json'), 'feature': record(W / 'producer-feature-acceptance.json'),
          'products': record(W / 'products-acceptance.json'), 'arithmetic': verify_arithmetic(),
          'producer_output_seals': {label: verify_output_seal(label) for label in P['producer_labels']},
          'executable': record(T / 'controls/original_candidate_controls'),
          'scope': 'Original eight controls only; all original filters, assertions, limits and provider/cache policy retained.'})
    print('Original eight controls admitted behind terminal compiler/profile/producer acceptance.', flush=True)

def review():
    verify_static(); verify_arithmetic(); verify_host_reuse()
    runtime = json.loads((W / 'runtime-context.json').read_text())
    assert runtime['producer_output_seals'] == {label: verify_output_seal(label) for label in P['producer_labels']}
    assert runtime['producer'] == record(W / 'producer-acceptance.json') and runtime['feature'] == record(W / 'producer-feature-acceptance.json')
    assert runtime['executable'] == record(T / 'controls/original_candidate_controls')
    actions = json.loads((W / 'actions.json').read_text())
    for index in range(1, 9):
        label = 'control-' + str(index); owner = accepted(label); join_runtime_inputs(label, owner)
        assert owner['limits'] == P['control_limits'][label]
        expected = set(P['provider_identities']); observed = set()
        births = {(row['pid'], row['saved_start_ticks']) for row in owner['observed_births']}
        for row in owner['provider_mappings']:
            assert (row['pid'], row['start_ticks']) in births and row['owned_session'] == owner['child_pid']
            for line in row['maps']:
                fields = line.split(None, 5); pin = P['provider_identities'][fields[5]]
                assert [int(v, 16) for v in fields[3].split(':')] == pin['device'] and int(fields[4]) == pin['inode']
                observed.add(fields[5])
        assert observed == expected
    for label in ('verify-controls', 'verify-endpoint-routes'): verify_output_seal(label)
    controls = json.loads((W / 'control-acceptance.json').read_text()); routes = json.loads((W / 'endpoint-route-acceptance.json').read_text())
    assert controls['passed'] and len(controls['original_controls']) == 8 and routes['passed']
    assert controls['timestamp_schema']['complete_readbacks'] == 32 and controls['timestamp_schema']['ordered_markers'] == 51
    assert controls['timestamp_schema']['actual_module_command_joins'] == 54
    assert routes['complete_readbacks'] == 56 and len(routes['routes']) == 4
    write(W / 'prerequisites-result.json', {'schema':'sirius-fp64-estimate-original-controls-result-v1','passed': True, 'original_eight_controls_passed': True,
          'profile': record(W / 'profile.json'), 'arithmetic': verify_arithmetic(),
          'control_acceptance': record(W / 'control-acceptance.json'), 'endpoint_route_acceptance': record(W / 'endpoint-route-acceptance.json'),
          'production_adopted': False, 'benefit_executed': False,
          'scope': 'Original primitive arithmetic and eight integrated prerequisites for explicit17-module exact-verified Float64 estimate diagnostic profile. Original20+2 adoption gates remain unexecuted; no full-science/frame/release qualification.'})
    print('Original eight controls and complete route/readback joins accepted; no benefit verdict.', flush=True)

if __name__ == '__main__':
    if sys.argv[1] == 'phase': phase(sys.argv[2])
    else: {'pin-context': pin_context, 'pin-runtime': pin_runtime, 'review': review}[sys.argv[1]]()
