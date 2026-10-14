"""Bounded original default-profile products, producer and eight-control owner actions."""
from pathlib import Path
import hashlib, json, re, shlex, shutil, subprocess, sys
sys.dont_write_bytecode = True
W = Path(__file__).resolve().parent
R = W.parents[1]
T = R / 'bin/linux-gcc/portable-multiply64-controls-o0'
P = json.loads((W / 'packet.json').read_text())

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

def verify_source_and_compiler():
    archive = Path(P['compiler_archive'])
    assert content(archive / 'manifest.json')['sha256'] == P['compiler_manifest_sha256']
    manifest = json.loads((archive / 'manifest.json').read_text())
    for name, pin in manifest['files'].items(): assert content(archive / name) == pin, name
    original = archive / 'files/out/portable-multiply64-compiler'
    accepted_folder(original / 'compile-endpoint-normal')
    result = json.loads((original / 'compiler-result.json').read_text())
    approval = json.loads((original / 'compiler-acceptance.json').read_text())
    assert approval['passed'] and approval['compiler_owner']['sha256'] == content(original / 'compile-endpoint-normal/owner.json')['sha256']
    assert approval['compiler_result']['sha256'] == content(original / 'compiler-result.json')['sha256']
    assert result['candidate_multiply_graph']['complete_discarded_bit_jam']
    assert result['original_divide_unpack_round_normalized_bodies_equal']
    prototype = archive / 'files/out/portable-multiply64-screen'
    assert (prototype / 'candidate/src/sirius/kernels/portable_binary32.h').read_bytes() == (W / 'candidate-header.h').read_bytes()
    assert (W / 'canonical-header.h').read_bytes() == (R / 'src/sirius/kernels/portable_binary32.h').read_bytes()
    for path in (W / 'stock/src/sirius/kernels').iterdir():
        assert path.read_bytes() == (R / 'src/sirius/kernels' / path.name).read_bytes()
        expected = W / 'candidate-header.h' if path.name == 'portable_binary32.h' else path
        assert (W / 'candidate/src/sirius/kernels' / path.name).read_bytes() == expected.read_bytes()
    return record(original / 'compiler-acceptance.json')

def verify_arithmetic():
    # No owner/context/compiler admission until the actual immutable numeric archive exists and is pinned.
    assert re.fullmatch(r'[0-9a-f]{64}', P['arithmetic_manifest_sha256']), 'original multiply64 arithmetic is still pending'
    assert re.fullmatch(r'[0-9a-f]{40}', P['arithmetic_evidence_commit'])
    archive = Path(P['arithmetic_archive'])
    assert content(archive / 'manifest.json')['sha256'] == P['arithmetic_manifest_sha256']
    manifest = json.loads((archive / 'manifest.json').read_text())
    for name, pin in manifest['files'].items(): assert content(archive / name) == pin, name
    work = archive / 'files/out' / P['arithmetic_work_name']
    result = json.loads((work / 'arithmetic-result.json').read_text())
    assert content(work / 'arithmetic-result.json') == P['arithmetic_result_pin']
    assert content(work / 'independent-terminal-review.json') == P['arithmetic_terminal_review_pin']
    assert content(archive / 'result.json') == P['arithmetic_aggregate_result_pin']
    terminal = json.loads((work / 'independent-terminal-review.json').read_text())
    assert terminal['passed'] and terminal['terminal_owners'] == 10 and terminal['output_seals'] == 8
    assert terminal['root_rechecked_pins']['arithmetic-result.json'] == P['arithmetic_result_pin']
    assert terminal['complete_candidate_output_words'] == 697072 and terminal['exact_raw_word_mismatches'] == 0
    assert terminal['full_candidate_original_expected_and_historical_stock_equal']
    assert terminal['host_total_original_records'] == 174268 and terminal['host_original_residual_records'] == 53524
    assert result['schema'] == 'sirius-multiply64-original-arithmetic-result-v1'
    assert result['passed'] and result['original_records'] == 174268
    assert result['original_operation_count'] == 14 and result['complete_output_words_per_provider'] == 697072
    assert result['host_explicit64_and_both_full_gpu_outputs_equal_original_expected']
    assert result['stock_result_reused'] and result['candidate_macro'] == 'SIRIUS_PB32_MULTIPLY64=1'
    assert result['original_corpus']['operation_cardinalities'] == P['arithmetic_expected_operation_cardinalities']
    candidate = result['modules']['candidate']; stock = result['modules']['stock']
    assert candidate['capabilities'] == ['Int64', 'Shader'] and candidate['integer_widths'] == ['32', '64']
    assert stock['capabilities'] == ['Shader'] and stock['integer_widths'] == ['32']
    assert candidate['divide_widths'] == stock['divide_widths'] == [32]
    assert candidate['divide_loops'] == stock['divide_loops'] == 1
    assert candidate['original_divide_unpack_round_normalized_bodies_equal']
    shape = candidate['multiply_graph']
    assert shape['upper_bit'] == 47 and shape['shift_values'] == [20, 21]
    assert shape['normalized_unsigned_operands'] and shape['complete_discarded_bit_jam'] and shape['same_common_round_call']
    assert result['stock_arithmetic_receipt']['all737_static_inputs_and17_terminal_owners_joined']
    assert result['candidate_shader'] == candidate['module']
    actual_module = archive / 'files' / Path(result['candidate_shader']['resolved_path']).relative_to(R)
    assert content(actual_module) == {k: result['candidate_shader'][k] for k in ('bytes', 'sha256')}
    assert (work / 'candidate/src/sirius/kernels/portable_binary32.h').read_bytes() == (W / 'candidate-header.h').read_bytes()
    for label, pin in result['owners'].items():
        accepted_folder(work / label)
        assert content(work / label / 'owner.json') == {k: pin[k] for k in ('bytes', 'sha256')}
    original_stock = Path(P['stock_arithmetic_archive']) / 'files/out/portable-divide64-arithmetic'
    assert content(Path(P['stock_arithmetic_archive']) / 'manifest.json')['sha256'] == P['stock_arithmetic_manifest_sha256']
    expected = (work / 'expected.bin').read_bytes()
    assert expected == (original_stock / 'expected.bin').read_bytes() == (original_stock / 'execute-gpu-stock/actual.bin').read_bytes()
    assert (work / 'execute-gpu-candidate/actual.bin').read_bytes() == expected
    return record(work / 'arithmetic-result.json')

def pin_context():
    verify_static(); verify_source_and_compiler(); verify_host_reuse(); verify_arithmetic()
    verify_completed_producers()
    assert not T.exists() and not (W / 'operator-context.json').exists()
    assert not subprocess.check_output(['git', 'status', '--porcelain=v1'], cwd=R)
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=R, text=True).strip()
    assert revision == P['execution_revision']
    assert (W / 'actions.json').read_bytes() == (W / 'actions.unpinned.json').read_bytes()
    actions = json.loads((W / 'actions.unpinned.json').read_text())
    for action in actions.values():
        assert action['revision'] == 'PENDING_OPERATOR_CLEAN_CONTEXT'
        action['revision'] = revision
    (W / 'actions.json').write_text(json.dumps(actions, sort_keys=True, indent=2) + '\n')
    write(W / 'operator-context.json', {'revision': revision, 'actions': record(W / 'actions.json'),
          'source_freeze': record(W / 'source-freeze.json'), 'preparation': record(W / 'preparation.json')})
    print('Pinned original optional64 profile; no compiler/provider launched.', flush=True)

def verify_completed_producers():
    archive = Path(P['completed_archive']); history = Path(P['completed_work'])
    assert content(archive / 'manifest.json') == P['completed_manifest_pin']
    manifest = json.loads((archive / 'manifest.json').read_text())
    for name, pin in manifest['files'].items(): assert content(archive / name) == pin, name
    old_packet = json.loads(Path(P['completed_original_packet']).read_text())
    old_actions = json.loads(Path(P['completed_original_actions']).read_text())
    frozen = json.loads((history / 'source-freeze.json').read_text())['pins']
    assert set(frozen) <= set(P['completed_input_joins'])
    for original, pin in frozen.items(): assert P['completed_input_joins'][original]['pin'] == pin
    receipts = {}
    for label in P['completed_labels']:
        owner = accepted_folder(history / label)
        before = json.loads((history / label / 'inputs-before.json').read_text())
        assert owner['argv'] == old_actions[label]['argv'] and owner['limits'] == old_actions[label]['limits']
        for key, pin in before.items():
            row = P['completed_input_joins'][str(R / key)]
            assert row['pin'] == pin
            actual = Path(row['actual'])
            assert content(actual) == {k: pin[k] for k in ('bytes', 'sha256')}
            if not actual.is_relative_to(archive): assert record(actual) == pin
        seal = json.loads((history / label / 'output-seal.json').read_text())
        assert seal['label'] == label and seal['recipe'] == old_packet['recipes'][label]
        for field, path in [('source_freeze', history / 'source-freeze.json'),
                            ('operator_context', history / 'operator-context.json'),
                            ('inputs_before', history / label / 'inputs-before.json')]:
            assert content(path) == {k: seal[field][k] for k in ('bytes', 'sha256')}
        assert set(seal['outputs']) == set(old_packet['outputs'][label])
        for original, pin in seal['outputs'].items():
            source = archive / 'files' / Path(original).relative_to(R)
            assert content(source) == {k: pin[k] for k in ('bytes', 'sha256')}
        receipts[label] = {'owner': record(history / label / 'owner.json'), 'output_seal': record(history / label / 'output-seal.json')}
    failure = P['original_product_reader_failure']
    stopped = json.loads(Path(failure['owner']).read_text())
    assert stopped['status'] == 'terminal' and stopped['passed'] is False and stopped['child_exit'] == 1
    assert stopped['stop_reason'] is None and stopped['cleanup_errors'] == [] and stopped['remaining_owned_processes'] == []
    assert stopped['source_unchanged'] and stopped['whole_inputs_unchanged'] and all(row['absent'] for row in stopped['observed_births'])
    assert content(failure['stderr']) == failure['stderr_pin']
    for directory in ('stock', 'candidate'):
        for original in (history / directory).rglob('*'):
            if original.is_file(): assert original.read_bytes() == (W / directory / original.relative_to(history / directory)).read_bytes()
    return {'passed': True, 'archive_manifest': record(archive / 'manifest.json'), 'completed_three_owners': receipts,
            'old_product_reader_owner_preserved_false': record(failure['owner']), 'source_input_joins': len(P['completed_input_joins']),
            'scope': 'Completed stock39/candidate39 compiler products and canonical selfcheck reused exactly. No compiler, validator or policy test rerun. Stopped product reader remains failed; new full product reader is independent.'}

def reuse_completed_products():
    proof = verify_completed_producers(); archive = Path(P['completed_archive'])
    for mode in ('stock', 'candidate'):
        source = archive / 'files/bin/linux-gcc/portable-multiply64-controls' / mode
        target = T / mode
        assert not target.exists(); target.mkdir(parents=True)
        names = sorted(path.name for path in source.iterdir()); assert len(names) == 39
        for name in names:
            shutil.copyfile(source / name, target / name)
            assert (source / name).read_bytes() == (target / name).read_bytes()
    proof['copied_files'] = 78
    write(W / 'completed-products-reuse.json', proof)


def products_review():
    from elf_arrays import arrays
    from interface_signature import interface_signature
    from multiply_profile import multiply_shape, normalized_function
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
            caps_stock = re.findall(r'OpCapability (\S+)', stock)
            caps = re.findall(r'OpCapability (\S+)', candidate)
            assert 'Int64' not in caps_stock and sorted(caps) == sorted(caps_stock + ['Int64']), name
            assert re.findall(r'OpTypeInt 64 ([01])', candidate) == ['0'], name
            assert re.findall(r'OpTypeFloat (\d+)', candidate) == re.findall(r'OpTypeFloat (\d+)', stock), name
            assert not re.search(r'\b(?:RelaxedPrecision|FPFastMathMode|FPFastMathDefault|OpFmaKHR)\b', candidate)
            camera = name.startswith(('retained_camera_portable', 'retained_ray_camera_portable'))
            native_o0 = name.endswith('_fp64.spv') and '_portable' not in name
            footprint = multiply_shape(candidate, camera=camera, allow_o0=native_o0)
            for function in ('PB32Unpack', 'PB32Round', 'PB32Divide', 'PB32WideProduct', 'PB32ProductResidual',
                             'OrderedAdd32', 'OrderedSubtract32', 'OrderedMultiply32', 'OrderedDivide32'):
                present = [bool(re.search(r'OpName %\S+ "' + function + '"', text)) for text in (stock, candidate)]
                assert present[0] == present[1], (name, function)
                if present[0]: assert normalized_function(stock, function) == normalized_function(candidate, function), (name, function)
            assert re.findall(r'= (OpF(?!unction)\S+) ', stock) == re.findall(r'= (OpF(?!unction)\S+) ', candidate), name
            wide = re.search(r'(%\S+) = OpTypeInt 64 0', candidate)[1]
            assert not re.search(r'Op(?:U|S)(?:Div|Mod|Rem) ' + re.escape(wide) + r' ', candidate)
        rows.append({'module': name, 'changed': changed, 'stock': record(paths['stock']), 'candidate': record(paths['candidate']),
                     'complete_public_shared_type_layout_modes_and_ordered_barriers_equal': True, 'multiply_footprint': footprint})
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
    return {'passed': True, 'completed_products_reuse': record(W / 'completed-products-reuse.json'), 'reader_controls': record(W / 'parser-controls-result.json'), 'modules': rows, 'headers': headers, 'changed_modules': 24, 'unchanged_modules': 11,
            'main_arrays_changed': 24, 'main_arrays_unchanged': 16, 'owner_arrays_unchanged': 2,
            'all_seven_tables_and_nonarray_ABI_constants_equal': True,
            'scope': 'Complete compiler/profile/ABI evidence only; no integrated numeric or benefit verdict.'}

def prepare_producer():
    assert json.loads((W / 'products-acceptance.json').read_text())['passed']
    verify_source_and_compiler()
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
    assert producer['unchanged_array_count'] == 16
    for name in ('kTransportShader', 'kEndpointShader'):
        assert producer['actual_ELF_owner_arrays']['baseline'][name]['sha256'] == producer['actual_ELF_owner_arrays']['candidate'][name]['sha256']
    return {'passed': True, 'all28device_definitions_in_object_and_actual_candidate': True,
            'old_vulkan_and_compute_archive_members_absent': True, 'link_uses_default_archive_extraction': True,
            'feature_host': verify_host_reuse(), 'nm_outputs': symbols,
            'scope': 'Candidate linked to identical admitted diagnostic device object. Prior public feature arguments remain historical; original baseline ELF is unchanged and is not claimed to enable Int64.'}

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
    if label == 'reuse-completed-products': reuse_completed_products()
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
    write(W / 'prerequisites-result.json', {'passed': True, 'original_eight_controls_passed': True,
          'profile': record(W / 'profile.json'), 'arithmetic': verify_arithmetic(),
          'control_acceptance': record(W / 'control-acceptance.json'), 'endpoint_route_acceptance': record(W / 'endpoint-route-acceptance.json'),
          'production_adopted': False, 'benefit_executed': False,
          'scope': 'Original primitive arithmetic and eight integrated prerequisites for explicit24-module diagnostic profile. Original20+2 adoption gates remain unexecuted; no full-science/frame/release qualification.'})
    print('Original eight controls and complete route/readback joins accepted; no benefit verdict.', flush=True)

if __name__ == '__main__':
    if sys.argv[1] == 'phase': phase(sys.argv[2])
    else: {'pin-context': pin_context, 'pin-runtime': pin_runtime, 'review': review}[sys.argv[1]]()
