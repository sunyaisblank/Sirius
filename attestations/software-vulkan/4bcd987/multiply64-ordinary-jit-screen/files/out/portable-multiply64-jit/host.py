"""Reuse accepted stock/Int64 host and admit one candidate-only ordinary capture."""
from pathlib import Path
import hashlib, importlib.util, json, os, re, shlex, shutil, subprocess, sys
sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[1]
BUILD = ROOT / 'bin/linux-gcc/portable-multiply64-jit'
P = json.loads((WORK / 'inputs.json').read_text())

def record(path):
    path = Path(path)
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': path.stat().st_size, 'sha256': digest, 'resolved_path': str(path.resolve(strict=True))}

def content(path):
    return {key: value for key, value in record(path).items() if key != 'resolved_path'}

def load(path):
    return json.loads(Path(path).read_text())

def write(path, value):
    path = Path(path)
    assert not path.exists(), path
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n')

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result

def canonical_json(value):
    # Exact durable JSON comparison from the accepted separate pair reader.
    # Normalize integer object keys only through their actual JSON persistence.
    normalized = json.loads(json.dumps(value, allow_nan=False))
    return json.dumps(normalized, sort_keys=True, separators=(',', ':'), allow_nan=False)

def verify_static():
    for name, pin in load(WORK / 'source-freeze.json')['pins'].items():
        assert record(name) == pin, name

def manifest(archive, expected):
    archive = Path(archive)
    assert content(archive / 'manifest.json')['sha256'] == expected
    result = load(archive / 'manifest.json')
    for name, pin in result['files'].items():
        path = archive / name
        assert not path.is_symlink() and content(path) == pin, path
    return result

def joined(original, pin, archive):
    original = Path(original)
    if not original.is_absolute(): original = ROOT / original
    relative = original.relative_to(ROOT) if original.is_relative_to(ROOT) else None
    saved = Path(archive) / 'files' / relative if relative is not None else None
    actual = saved if saved is not None and saved.is_file() else original
    assert content(actual) == {key: pin[key] for key in ('bytes', 'sha256')}, original
    if actual == original and 'resolved_path' in pin:
        assert str(actual.resolve(strict=True)) == pin['resolved_path']
    return actual

def accepted(folder):
    folder = Path(folder)
    owner = load(folder / 'owner.json')
    assert owner['status'] == 'terminal' and owner['passed'] is True and owner['child_exit'] == 0
    assert owner['stop_reason'] is None and owner['cleanup_errors'] == []
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['remaining_owned_processes'] == [] and owner['observed_births']
    assert all(row['absent'] for row in owner['observed_births'])
    assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
    return owner

def verify_host_reuse():
    archive = Path(P['historical_archive'])
    manifest(archive, P['historical_manifest_sha256'])
    history = Path(P['historical_work'])
    frozen = load(history / 'source-freeze.json')['pins']
    assert len(frozen) == len(P['historical_frozen_joins']) == 605
    for row in P['historical_frozen_joins']:
        assert row['pin'] == frozen[row['original']]
        assert joined(row['original'], row['pin'], archive) == Path(row['actual'])
    for phase in ('compile-device', 'compile-probe', 'link'):
        accepted(history / phase)
    for phase, stem in (('compile-device', 'int64-vulkan_device'), ('compile-probe', 'short_input_jit_probe')):
        seal = load(history / phase / 'output-seal.json')
        assert seal['phase'] == phase
        for key, relative in (('source_freeze', 'source-freeze.json'), ('operator_context', 'operator-context.json'),
                              ('inputs_before', phase + '/inputs-before.json')):
            assert content(history / relative) == {field: seal[key][field] for field in ('bytes', 'sha256')}
        for name, pin in seal['outputs'].items(): joined(name, pin, archive)
        before = load(history / phase / 'inputs-before.json')
        normalized = {str((ROOT / key).resolve()): value for key, value in before.items()}
        dependency = archive / 'files/bin/linux-gcc/portable-divide64-jit' / (stem + '.d')
        for name in shlex.split(dependency.read_text().replace('\\\n', '').split(':', 1)[1]):
            original = str((ROOT / name).resolve())
            assert original in normalized
            joined(original, normalized[original], archive)
    link = load(history / 'link-review.json')
    assert link['passed'] and link['all_required_definitions_in_private_object_and_executable']
    assert link['private_object_before_original_archive'] and not link['original_vulkan_archive_member_extracted']
    assert len(link['required_strong_sirius_definitions']) == 28
    assert joined(link['executable']['resolved_path'], link['executable'], archive) == Path(P['canary'])
    joined(link['link_map']['resolved_path'], link['link_map'], archive)
    for phase, pin in link['compile_output_seals'].items():
        assert content(history / phase / 'output-seal.json') == {key: pin[key] for key in ('bytes', 'sha256')}
    before = load(history / 'link/inputs-before.json')
    for arg in load(history / 'host-recipes.json')['link']:
        if Path(arg).is_absolute():
            path = Path(arg)
            key = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
            if key in before: joined(path, before[key], archive)
    stock = history / 'capture-stock-endpoint-normal'
    accepted(stock)
    stock_code = load(stock / 'actual-code-review.json')
    stock_features = load(stock / 'feature-review.json')
    assert stock_code['actual_loaded_module']['sha256'] == '981bc85d13e168c6adf72c99adfe169157593dda1bee46e78ff867fc6b199ede'
    assert stock_features['actual_public_create_arguments_enable_shaderInt64']
    assert stock_features['successful_logical_device_used_for_actual_complete_shader']
    assert len(stock_features['features']) == 55 and stock_features['features']['shaderInt64'] == 1
    return {'historical_manifest': record(archive / 'manifest.json'), 'original_compile_device': record(history / 'compile-device/owner.json'),
            'original_compile_probe': record(history / 'compile-probe/owner.json'), 'original_link': record(history / 'link/owner.json'),
            'link_review': record(history / 'link-review.json'), 'canary': record(P['canary']),
            'historical_stock_owner': record(stock / 'owner.json'), 'historical_stock_code': record(stock / 'actual-code-review.json'),
            'historical_stock_features': record(stock / 'feature-review.json'),
            'scope': 'Exact accepted source/object/link/tool/library/host and complete stock capture reuse. Old public feature arguments remain historical; no stock recapture or recompile.'}

def verify_compiler():
    archive = Path(P['compiler_archive'])
    manifest(archive, P['compiler_manifest_sha256'])
    work = Path(P['compiler_work'])
    accepted(work / 'compile-endpoint-normal')
    for name, expected in (('compiler-acceptance.json', 'compiler_acceptance_pin'),
                            ('compiler-result.json', 'compiler_result_pin'),
                            ('compile-endpoint-normal/owner.json', 'compiler_owner_pin')):
        assert content(work / name) == P[expected]
    result = load(work / 'compiler-result.json'); approval = load(work / 'compiler-acceptance.json')
    assert approval['passed'] and approval['compiler_owner']['sha256'] == P['compiler_owner_pin']['sha256']
    assert approval['compiler_result']['sha256'] == P['compiler_result_pin']['sha256']
    assert result['candidate_multiply_graph']['complete_discarded_bit_jam']
    assert result['original_divide_unpack_round_normalized_bodies_equal']
    assert result['complete_descriptor_and_workgroup_type_decorations_exact'] and result['execution_modes_exact']
    assert result['ordered_barrier_operands_exact'] and not result['runtime_or_numerical_admission']
    for values in result['modules'].values():
        for pin in (values['module'], values['assembly'], *values['bodies'].values()): joined(pin['resolved_path'], pin, archive)
    for key in ('private_inputs', 'external_inputs'):
        for name, pin in load(work / 'source-freeze.json')[key].items():
            original = ROOT / 'out/portable-multiply64-compiler' / name if key == 'private_inputs' else Path(name)
            joined(original, pin, archive)
    action = load(work / 'actions.json')['compile-endpoint-normal']
    before = load(work / 'compile-endpoint-normal/inputs-before.json')
    for name in action['frozen_inputs']:
        path = Path(name); key = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
        joined(path, before[key], archive)
    assert joined(P['candidate_module_original_path'], result['modules']['candidate']['module'], archive) == Path(P['candidate_module'])
    assert content(P['candidate_module']) == P['candidate_module_pin']
    return {'archive_manifest': record(archive / 'manifest.json'), 'owner': record(work / 'compile-endpoint-normal/owner.json'),
            'acceptance': record(work / 'compiler-acceptance.json'), 'result_with_all_module_body_output_pins': record(work / 'compiler-result.json'),
            'candidate_module': record(P['candidate_module']),
            'scope': 'Immutable complete compiler result output pins, acceptance, terminal owner, source and consumed-input joins. No separate output-seal.json existed or is fabricated.'}

def verify_providers():
    for name, pin in P['providers'].items():
        path = Path(name); current = path.stat()
        assert record(path) == {key: pin[key] for key in ('bytes', 'sha256', 'resolved_path')}
        assert [os.major(current.st_dev), os.minor(current.st_dev)] == pin['device'] and current.st_ino == pin['inode']

def pin_runtime():
    verify_static(); host = verify_host_reuse(); compiler = verify_compiler(); verify_providers()
    assert not BUILD.exists() and not (WORK / 'module-inputs').exists()
    assert not any((WORK / name).exists() for name in ('runtime-context.json', 'capture-actions.json', 'reuse-acceptance.json', 'capture-candidate-endpoint-normal'))
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip() == P['execution_revision']
    assert not subprocess.check_output(['git', 'status', '--porcelain=v1'], cwd=ROOT)
    BUILD.mkdir()
    executable = BUILD / 'short_input_jit_probe'
    shutil.copyfile(P['canary'], executable); executable.chmod(0o755)
    assert content(executable) == content(P['canary'])
    (WORK / 'module-inputs').mkdir()
    candidate = WORK / 'module-inputs/candidate-endpoint-normal.spv'
    shutil.copyfile(P['candidate_module'], candidate)
    assert content(candidate) == P['candidate_module_pin']
    actions = load(WORK / 'capture-actions.unpinned.json')
    assert list(actions) == ['capture-candidate-endpoint-normal']
    action = actions['capture-candidate-endpoint-normal']
    assert action['revision'] == 'PENDING_OPERATOR_CLEAN_CONTEXT'
    action['revision'] = P['execution_revision']
    write(WORK / 'capture-actions.json', actions)
    write(WORK / 'reuse-acceptance.json', {'host': host, 'compiler': compiler,
          'stock_recaptured': False, 'host_recompiled': False, 'passed': True})
    write(WORK / 'runtime-context.json', {'revision': P['execution_revision'], 'host': record(executable),
          'modules': {'candidate': record(candidate)}, 'capture_actions': record(WORK / 'capture-actions.json'),
          'source_freeze': record(WORK / 'source-freeze.json'), 'reuse_acceptance': record(WORK / 'reuse-acceptance.json'),
          'scope': 'One candidate-only original short-input canary. Complete stock/Int64 host evidence reused; no arithmetic, performance, cold-frame or qualification claim.'})
    print('Candidate-only runtime inputs pinned; no provider process launched.', flush=True)

def verify_new_maps(folder, owner, code):
    verify_providers()
    pid = load(folder / 'inferior.json')['pid']
    births = {row['pid']: row['saved_start_ticks'] for row in owner['observed_births']}
    assert pid in births
    expected = set(P['providers'])
    def check(lines):
        found = set()
        for line in lines:
            fields = line.split(None, 5)
            if len(fields) != 6 or fields[5] not in expected: continue
            pin = P['providers'][fields[5]]
            assert [int(value, 16) for value in fields[3].split(':')] == pin['device']
            assert int(fields[4]) == pin['inode']; found.add(fields[5])
        return found
    observations = [row for row in owner['provider_mappings'] if row['pid'] == pid]
    assert observations and all(row['start_ticks'] == births[pid] and row['owned_session'] == owner['child_pid'] for row in observations)
    assert set().union(*(check(row['maps']) for row in observations)) == expected
    bindings = load(folder / 'bindings.json')
    events = [row for row in bindings['jit_events'] if row.get('file') == code['actual_entry']['file']]
    assert len(events) == 1 and check(events[0]['maps_at_first_call'].splitlines()) == expected

def review():
    verify_static(); verify_host_reuse(); verify_compiler()
    folder = WORK / 'capture-candidate-endpoint-normal'
    owner = accepted(folder); context = load(WORK / 'runtime-context.json')
    action = load(WORK / 'capture-actions.json')['capture-candidate-endpoint-normal']
    assert owner['source'] == {'revision': P['execution_revision'], 'status': ''}
    assert owner['argv'] == [arg.replace('{folder}', str(folder)) for arg in action['argv']]
    assert owner['limits'] == action['limits']
    assert context['capture_actions'] == record(WORK / 'capture-actions.json')
    assert context['source_freeze'] == record(WORK / 'source-freeze.json')
    assert context['host'] == record(BUILD / 'short_input_jit_probe')
    assert context['modules']['candidate'] == record(WORK / 'module-inputs/candidate-endpoint-normal.spv')
    before = load(folder / 'inputs-before.json')
    for name in action['frozen_inputs']:
        path = Path(name); key = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
        assert before[key] == record(path), name
    analyzer = module('original_jit_analysis', WORK / 'analyze_capture.py')
    feature_reader = module('original_feature_reader', WORK / 'read_features.py')
    code = analyzer.inspect(folder); features = feature_reader.inspect(folder)
    stock_folder = Path(P['historical_work']) / 'capture-stock-endpoint-normal'
    stock_code = load(stock_folder / 'actual-code-review.json')
    assert canonical_json(analyzer.inspect(stock_folder)) == canonical_json(stock_code)
    stock_features = load(stock_folder / 'feature-review.json')
    assert feature_reader.inspect(stock_folder) == stock_features
    layout = load(WORK / 'feature-layout.json')
    assert len(layout['features_fields']) == len(set(layout['features_fields'])) == 55
    assert set(features['features']) == set(layout['features_fields']) and features == stock_features
    assert code['actual_loaded_module']['sha256'] == P['candidate_module_pin']['sha256']
    assert code['actual_loaded_module']['bytes'] == P['candidate_module_pin']['bytes']
    assert content(BUILD / 'short_input_jit_probe') == content(P['canary'])
    verify_new_maps(folder, owner, code)
    write(folder / 'actual-code-review.json', code)
    write(folder / 'feature-review.json', features)
    write(WORK / 'capture-acceptance.json', {'passed': True, 'candidate_owner': record(folder / 'owner.json'),
          'candidate_code_review': record(folder / 'actual-code-review.json'), 'candidate_feature_review': record(folder / 'feature-review.json'),
          'historical_stock_code_review': record(stock_folder / 'actual-code-review.json'),
          'historical_stock_feature_review': record(stock_folder / 'feature-review.json'),
          'complete_durable_stock_report_reproduced': True, 'complete_55_feature_queue_extension_order_and_chain_equal': True,
          'fresh_inferior_and_called_JIT_maps_join_all_three_provider_files': True,
          'scope': 'One new candidate loaded/called ordinary IR/JIT capture with exact reused host and historical stock. New candidate feature arguments are fresh; stock feature arguments are historical. No retained multiplication, arithmetic equivalence, cost benefit or qualification verdict; lowering review remains separate.'})
    print('Candidate owner, ordinary code/relocations/features and fresh maps accepted; cost premise unjudged.', flush=True)

if __name__ == '__main__':
    {'pin-runtime': pin_runtime, 'review': review}[sys.argv[1]]()
