"""Restore the accepted guarded-FP64-estimate producer pair for the unchanged original trial."""
from pathlib import Path
import ast, hashlib, json, re, shutil, subprocess
import xml.etree.ElementTree as ET

R = Path.cwd().resolve()
Q = Path(__file__).resolve().parent
W = Q / 'benefit'
T = R / 'bin/linux-gcc/portable-divide-fp64-estimate-benefit'
for name, pin in json.loads((Q / 'source-freeze.json').read_text())['pins'].items():
    path = Path(name)
    data = path.read_bytes()
    assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == pin, name
P = json.loads((Q / 'prerequisite-pin.json').read_text())
# Fail before any trial output or owner admission while preservation is pending.
assert P['admitted'] is True, 'immutable eight-control PASS and final source admission pending'
assert re.fullmatch('[0-9a-f]{64}', P['manifest_sha256'])
assert re.fullmatch('[0-9a-f]{40}', P['evidence_commit'])
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == P['execution_revision']
assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
subprocess.run(['git', 'cat-file', '-e', P['evidence_commit'] + '^{commit}'], check=True)
assert not T.exists() and not W.exists() and not (Q / 'preparation.json').exists()

def identity(path):
    path = Path(path)
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': path.stat().st_size, 'sha256': digest}

def write(path, value):
    assert not path.exists(), path
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')

archive = R / P['archive']
manifest_path = archive / 'manifest.json'
assert identity(manifest_path)['sha256'] == P['manifest_sha256']
manifest = json.loads(manifest_path.read_text())
old_work = 'out/portable-divide-fp64-estimate-controls-continuation/'
old_bin = 'bin/linux-gcc/portable-divide-fp64-estimate-controls-continuation/'
old_stage = old_bin + 'controls/'
old_roots = [(R / name).resolve() for name in (old_work, old_bin)]
archive_inputs = {str(manifest_path): identity(manifest_path)}
current_inputs = {}

def preserved(relative):
    assert not Path(relative).is_absolute() and '..' not in Path(relative).parts
    path = archive / 'files' / relative
    assert not path.is_symlink(), path
    assert identity(path) == manifest['files']['files/' + relative], relative
    archive_inputs[str(path)] = identity(path)
    return path

def joined(original, pin):
    original = Path(original)
    if not original.is_absolute():
        original = R / original
    relative = original.relative_to(R).as_posix() if original.is_relative_to(R) else None
    if relative and relative.startswith((old_work, old_bin)):
        actual = preserved(relative)
        assert pin.get('resolved_path', str(original)) == str(original), original
    else:
        actual = original
        assert pin.get('resolved_path', str(actual.resolve(strict=True))) == str(actual.resolve(strict=True))
        current_inputs[str(actual)] = identity(actual)
    assert identity(actual) == {key: pin[key] for key in ('bytes', 'sha256')}, original
    return actual

packet = json.loads(preserved(old_work + 'packet.json').read_text())
actions = json.loads(preserved(old_work + 'actions.json').read_text())
context = json.loads(preserved(old_work + 'operator-context.json').read_text())
assert context['revision'] == P['historical_execution_revision']
assert identity(preserved(old_work + 'source-freeze.json'))['sha256'] == P['source_freeze_sha256']
assert identity(preserved(old_work + 'gate.py'))['sha256'] == P['gate_sha256']
for name, pin in json.loads(preserved(old_work + 'source-freeze.json').read_text())['pins'].items():
    joined(name, pin)
assert packet['changed_arrays'] == P['changed_arrays'] and len(P['changed_arrays']) == 17
assert packet['candidate_macro'] == P['candidate_macro'] == 'SIRIUS_PB32_DIVIDE_FP64_ESTIMATE=1'
assert packet['arithmetic_commit'] == P['arithmetic_commit']
assert packet['arithmetic_manifest_sha256'] == P['arithmetic_manifest_sha256']
assert packet['arithmetic_inner_result_sha256'] == P['arithmetic_inner_result_sha256']
assert packet['arithmetic_terminal_review_sha256'] == P['arithmetic_terminal_review_sha256']
for key in ('compiler_stop_commit', 'compiler_stop_manifest_sha256', 'profile_archive_commit',
            'profile_archive_manifest_sha256', 'profile_acceptance_sha256', 'profile_terminal_review_sha256'):
    assert packet[key] == P[key]

owners = {}
def accepted(label):
    folder = old_work + label + '/'
    owner_path = preserved(folder + 'owner.json')
    owner = json.loads(owner_path.read_text())
    assert owner['status'] == 'terminal' and owner['passed'] is True and owner['child_exit'] == 0, label
    assert owner['source'] == {'revision': context['revision'], 'status': ''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and owner['cleanup_errors'] == []
    assert owner['remaining_owned_processes'] == []
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    before_path = preserved(folder + 'inputs-before.json')
    assert before_path.read_bytes() == preserved(folder + 'inputs-after.json').read_bytes()
    for name, pin in json.loads(before_path.read_text()).items():
        joined(name, pin)
    expected = [arg.replace('{folder}', str(R / old_work / label)) for arg in actions[label]['argv']]
    assert owner['argv'] == expected and owner['limits'] == actions[label]['limits']
    owners[label] = identity(owner_path)
    return owner

producer_labels = packet['producer_labels']
labels = producer_labels + ['control-' + str(i) for i in range(1, 9)] + ['verify-controls', 'verify-endpoint-routes']
assert len(producer_labels) == 5 and len(labels) == len(set(labels)) == 15
for label in labels:
    accepted(label)
for label in producer_labels + ['verify-controls', 'verify-endpoint-routes']:
    seal = json.loads(preserved(old_work + label + '/output-seal.json').read_text())
    assert seal['label'] == label and seal['recipe'] == packet['recipes'][label]
    joined(R / old_work / 'source-freeze.json', seal['source_freeze'])
    joined(R / old_work / 'operator-context.json', seal['operator_context'])
    joined(R / old_work / label / 'inputs-before.json', seal['inputs_before'])
    assert set(seal['outputs']) == set(packet['outputs'][label])
    for name, pin in seal['outputs'].items():
        joined(name, pin)

result_path = preserved(old_work + 'prerequisites-result.json')
assert identity(result_path)['sha256'] == P['prerequisites_result_sha256']
result = json.loads(result_path.read_text())
assert result['schema'] == 'sirius-fp64-estimate-controls-continuation-result-v1'
assert result['fresh_owner_count'] == 15 and result['fresh_producer_count'] == 5
terminal_path = preserved(old_work + 'independent-terminal-review.json')
assert identity(terminal_path)['sha256'] == P['terminal_prerequisite_review_sha256']
assert json.loads(terminal_path.read_text())['passed'] is True
assert result['passed'] and result['original_eight_controls_passed']
assert result['production_adopted'] is False and result['benefit_executed'] is False
for key in ('profile', 'control_acceptance', 'endpoint_route_acceptance'):
    joined(result[key]['resolved_path'], result[key])
for key in ('result', 'terminal_review', 'manifest'):
    joined(result['arithmetic'][key]['resolved_path'], result['arithmetic'][key])
assert result['arithmetic']['all_original_static_and10owner_input_joins'] is True
# Reuse the accepted immutable receipt reader, without any new archive adapter.
receipt_path = preserved(old_work + 'prerequisite_receipts.py')
assert identity(receipt_path)['sha256'] == P['prerequisite_reader_sha256']
receipt_tree = ast.parse(receipt_path.read_text())
receipt_definitions = ast.Module(body=[n for n in receipt_tree.body if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef))], type_ignores=[])
receipt_reader = {}
exec(compile(receipt_definitions, 'sealed_accepted_prerequisite_reader', 'exec'), receipt_reader)
inherited = receipt_reader['verify'](packet, R)
assert inherited == result['inherited_prerequisites']
assert inherited['arithmetic'] == result['arithmetic']
assert len(inherited['historical_three_producers']) == 3
assert inherited['original_failed_wrapper_passed'] is False and inherited['original_failed_wrapper_child_exit'] == 1
controls = json.loads(preserved(old_work + 'control-acceptance.json').read_text())
routes = json.loads(preserved(old_work + 'endpoint-route-acceptance.json').read_text())
assert controls['passed'] and len(controls['original_controls']) == 8 and routes['passed']
schema = controls['timestamp_schema']
assert (schema['complete_readbacks'], schema['ordered_markers'], schema['actual_module_command_joins']) == (32, 51, 54)
assert schema['all_readbacks_exact_to_repaired_original_baseline']
assert routes['complete_readbacks'] == 56 and routes['complete_readbacks_exact_to_prior_original_control']
assert len(routes['routes']) == 4
assert {(row['declined_RTE32'], row['wide']) for row in routes['routes']} == {(d, w) for d in (False, True) for w in (False, True)}
for i in range(1, 9):
    label = 'control-' + str(i)
    owner = json.loads(preserved(old_work + label + '/owner.json').read_text())
    assert owner['provider_mappings']
    assert owner['limits'] == packet['control_limits'][label]
    row = controls['original_controls'][i - 1]
    assert row['label'] == label and row['owner'] == identity(preserved(old_work + label + '/owner.json'))
    xml_path = preserved(old_work + label + '/gtest.xml')
    assert row['xml'] == identity(xml_path)
    xml = ET.parse(xml_path).getroot()
    assert int(xml.get('tests', '-1')) == 1
    assert not any(xml.findall('.//' + tag) for tag in ('failure', 'error', 'skipped'))
    for node in [xml, *xml.findall('.//testsuite')]:
        assert all(int(node.get(key, '0')) == 0 for key in ('failures', 'errors', 'skipped', 'disabled'))
    cases = xml.findall('.//testcase')
    assert len(cases) == 1 and cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
    assert owner['argv'][1] == '--gtest_filter=' + cases[0].get('classname') + '.' + cases[0].get('name')

protected = json.loads(preserved(old_work + 'protected-before.json').read_text())
assert len(protected) == controls['protected_payloads_unchanged'] == 114
for name, pin in protected.items():
    joined(name, pin)
profile_result_path = joined(inherited['profile_result']['resolved_path'], inherited['profile_result'])
full_profile_result = json.loads(profile_result_path.read_text())
assert full_profile_result['passed']
products = full_profile_result['profile']
assert products['passed'] and len(products['modules']) == 35
assert (products['changed_modules'], products['unchanged_modules']) == (17, 18)
assert (products['main_arrays_changed'], products['main_arrays_unchanged'], products['owner_arrays_unchanged']) == (17, 23, 2)
assert products['all_seven_tables_and_nonarray_ABI_constants_equal']
producer = json.loads(preserved(old_work + 'producer-acceptance.json').read_text())
assert producer['passed'] and sorted(producer['changed_arrays']) == P['changed_arrays'] and producer['unchanged_array_count'] == 23
feature = json.loads(preserved(old_work + 'producer-feature-acceptance.json').read_text())
assert feature['passed'] and feature['all28device_definitions_in_object_and_actual_candidate']
assert feature['link_uses_default_archive_extraction'] and feature['old_vulkan_and_compute_archive_members_absent']
for pin in feature['feature_host'].values():
    if isinstance(pin, dict) and 'resolved_path' in pin:
        joined(pin['resolved_path'], pin)
runtime = json.loads(preserved(old_work + 'runtime-context.json').read_text())
for key in ('profile', 'producer', 'feature', 'executable'):
    joined(runtime[key]['resolved_path'], runtime[key])
assert runtime['arithmetic'] == result['arithmetic']
assert runtime['inherited_full_profile'] == inherited
assert runtime['producer_output_seals'] == {label: dict(identity(preserved(old_work + label + '/output-seal.json')), resolved_path=str(R / old_work / label / 'output-seal.json')) for label in producer_labels}

arithmetic = R / P['arithmetic_archive']
assert identity(arithmetic / 'manifest.json')['sha256'] == P['arithmetic_manifest_sha256']
current_inputs[str(arithmetic / 'manifest.json')] = identity(arithmetic / 'manifest.json')
subprocess.run(['git', 'cat-file', '-e', P['arithmetic_commit'] + '^{commit}'], check=True)
arith_work = arithmetic / 'files/out/portable-divide-fp64-estimate-arithmetic'
arith_result = arith_work / 'arithmetic-result.json'
assert identity(arith_result) == {k: result['arithmetic']['result'][k] for k in ('bytes', 'sha256')}
assert identity(arith_result)['sha256'] == P['arithmetic_inner_result_sha256']
assert identity(arithmetic / 'manifest.json') == {k: result['arithmetic']['manifest'][k] for k in ('bytes', 'sha256')}
ar = json.loads(arith_result.read_text())
assert ar['schema'] == 'sirius-fp64-estimate-original-arithmetic-result-v1'
assert ar['passed'] and ar['original_records'] == 174268 and ar['original_operation_count'] == 14
assert ar['complete_output_words_per_provider'] == 697072
assert ar['fresh_candidate_host_and_full_gpu_plus_reused_stock_equal_original_expected'] and ar['stock_result_reused']
assert ar['candidate_macro'] == P['candidate_macro']
assert ar['original_corpus']['division_branch_coverage']['optional_finite_nonzero_branch_records'] == 12263
assert {k: ar['candidate_header'][k] for k in ('bytes', 'sha256')} == P['candidate_header']
arith_header = arith_work / 'candidate/src/sirius/kernels/portable_binary32.h'
assert identity(arith_header) == P['candidate_header']
arith_terminal = arith_work / 'independent-terminal-review.json'
assert identity(arith_terminal) == {k: result['arithmetic']['terminal_review'][k] for k in ('bytes', 'sha256')}
assert identity(arith_terminal)['sha256'] == P['arithmetic_terminal_review_sha256']
assert json.loads(arith_terminal.read_text())['passed'] is True
for path in (arith_result, arith_header, arith_terminal):
    current_inputs[str(path)] = identity(path)

# Reuse only the original ELF-reader definitions, without running its footer.
tree = ast.parse(preserved(old_work + 'verify_producer.py').read_text())
definitions = ast.Module(body=[n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef))], type_ignores=[])
reader = {}
exec(compile(definitions, 'sealed_original_ELF_reader', 'exec'), reader)
for mode in ('baseline', 'candidate'):
    directory = 'stock' if mode == 'baseline' else 'candidate'
    generated = Path(packet['candidate_products']).parent / directory
    headers = {}
    for name in ('retained_kernels.h', 'retained_owner_kernels.h', 'retained_owner_layout.h'):
        original_header = generated / name
        profile_pin = full_profile_result['completed_producers']['sealed78_products'][str(original_header)]
        saved_header = joined(original_header, profile_pin)
        actual_header = saved_header if mode == 'baseline' else preserved(old_stage + 'headers/' + name)
        assert identity(actual_header) == identity(saved_header)
        headers[name] = actual_header
    executable = preserved(old_stage + 'original_' + mode + '_controls')
    assert identity(executable) == producer[mode + '_executable']
    assert identity(executable)['sha256'] == P[mode + '_sha256']
    actual, _ = reader['elf_arrays'](executable, reader['arrays'](headers['retained_kernels.h']))
    assert actual == producer['actual_ELF_arrays'][mode]
assert producer['link_owner'] == owners['link-candidate']
assert producer['link_map'] == identity(preserved(old_stage + 'candidate.map'))
assert producer['old_archive_compute_member_excluded']
for label, stem in (('compile-compute-sealed', 'compute-sealed'), ('compile-observer', 'observer')):
    proof = producer['compile_provenance'][label]
    assert proof['owner'] == owners[label]
    assert proof['object'] == identity(preserved(old_stage + stem + '.o'))
    for name, pin in proof['actual_dependencies'].items():
        joined(name, pin)

# All acceptance above precedes creation of the restored execution roots.
T.mkdir()
W.mkdir()
for name in ('producer-acceptance.json', 'producer-feature-acceptance.json', 'control-acceptance.json',
             'endpoint-route-acceptance.json', 'audited-observer.cpp', 'protected-before.json'):
    shutil.copyfile(preserved(old_work + name), Q / name)
assert identity(Q / 'audited-observer.cpp')['sha256'] == '399cd52cc6d9fc89e4b97ff63b06a54b06b4f80128f5a99c2c6d1278c62f3dc3'
executables = {}
for mode in ('baseline', 'candidate'):
    source = preserved(old_stage + 'original_' + mode + '_controls')
    target = T / source.name
    shutil.copyfile(source, target)
    target.chmod(0o755)
    assert identity(target) == producer[mode + '_executable']
    executables[mode] = {'path': str(target.relative_to(R)), 'source_revision': P['execution_revision'],
                         **identity(target), 'observer_source_sha256': identity(Q / 'audited-observer.cpp')['sha256']}
resources = []
for name in manifest['files']:
    relative = name.removeprefix('files/')
    if relative.startswith(old_stage + 'resources/'):
        source = preserved(relative)
        target = T / 'resources' / relative.removeprefix(old_stage + 'resources/')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        assert identity(target) == identity(source)
        resources.append(str(target))
assert len(resources) == 23
apparatus = Q / 'original-apparatus'
for name, pin in P['apparatus'].items():
    assert identity(apparatus / name) == pin
for name in ('software_contracts.py', 'compare_software.py', 'run_sequence.py', 'expected-readback-layout.json', 'run_owned.py', 'owned_guard.py'):
    shutil.copyfile(apparatus / name, W / name)
wrapper = (apparatus / 'run_limb_benefit_slot.py').read_text()
old_allowlist = 'bin/linux-gcc/computed-carry-exit-benefit/original_'
assert wrapper.count(old_allowlist) == 1
(W / 'run_limb_benefit_slot.py').write_text(wrapper.replace(old_allowlist, 'bin/linux-gcc/portable-divide-fp64-estimate-benefit/original_'))
for mode, name in (('baseline', 'baseline.json'), ('candidate', 'build-and-arrays.json')):
    write(W / name, {'arrays': producer['actual_ELF_arrays'][mode]})
shared = [*archive_inputs, *current_inputs, *resources]
shared += [str(Q / name) for name in ('prepare.py', 'prove_reuse.py', 'finish_trial.py', 'prerequisite-pin.json',
           'source-freeze.json', 'preparation.json', 'producer-reuse.json', 'producer-acceptance.json', 'producer-feature-acceptance.json',
           'control-acceptance.json', 'endpoint-route-acceptance.json', 'audited-observer.cpp', 'protected-before.json')]
shared += [str(W / name) for name in ('software_contracts.py', 'compare_software.py', 'run_sequence.py',
           'run_limb_benefit_slot.py', 'baseline.json', 'build-and-arrays.json', 'expected-readback-layout.json', 'retention-gate-frozen.json')]
assert not any(Path(name).resolve().is_relative_to(root) for name in shared for root in old_roots)
old_actions = json.loads((apparatus / 'actions.json').read_text())
new_actions = {}
for label in ('a1', 'b1', 'b2', 'a2'):
    mode = 'baseline' if label.startswith('a') else 'candidate'
    new_actions[label] = {'revision': P['execution_revision'], 'kind': 'test', 'expected_exit': 0,
        'argv': ['python3', '-B', str(W / 'run_limb_benefit_slot.py'), executables[mode]['path'], old_actions[label]['argv'][-1]],
        'limits': old_actions[label]['limits'], 'software_provider': True,
        'frozen_inputs': list(dict.fromkeys(shared + [str(R / executables[mode]['path'])])),
        'scope': 'One original fixed A/B/B/A slot, three repeats/family, unchanged20+2 gates. Full owned wall retained; device fixture SetUp excluded only as original. No pooling/subtraction/retry or full qualification.'}
    assert new_actions[label]['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
write(W / 'actions.json', new_actions)
freeze = json.loads((apparatus / 'retention-gate-frozen.json').read_text())
freeze.update(candidate_revision=P['execution_revision'], baseline_revision=P['execution_revision'],
              executables=executables, actions_sha256=identity(W / 'actions.json')['sha256'],
              scope='One original fixed guarded-FP64-estimate A1/B1/B2/A2 trial. Original stock host and accepted diagnostic Float64 host with historical unused Int64; candidate shaders require no Int64. Accepted17/23 plus2 unchanged owner arrays. Original measured windows, forced cacheoff and20+2 gates unchanged.')
assert freeze['sequence'] == ['a1', 'b1', 'b2', 'a2'] and freeze['within_slot_repeats'] == 3
assert (freeze['complete_call_minimum_improvement_percent'], freeze['every_inclusive_stage_minimum_improvement_percent'], freeze['observed_instance_api_minimum_improvement_percent']) == (5, -5, 0)
assert freeze['all20_gates_required']
assert W.parents[2] == R
write(W / 'retention-gate-frozen.json', freeze)
write(Q / 'preparation.json', {'source_revision': P['execution_revision'],
    'prior_manifest_sha256': P['manifest_sha256'], 'prerequisite_evidence_commit': P['evidence_commit'],
    'producer_prerequisites': {'actual40array_ELF_proof': True, 'actual2owner_array_ELF_proof': True,
        'same_observer_CPP': True, 'original8controls_passed': True, 'fifteen_fresh_terminal_owners': owners,
        'inherited_prerequisites': inherited,
        'candidate_diagnostic_Float64_host_with_historical_unused_Int64': True,
        'candidate_shaders_require_no_Int64': True, 'original_stock_host_unchanged': True},
    'executables': executables, 'changed_arrays': P['changed_arrays'],
    'current_byte_exact_dependency_joins': current_inputs, 'preserved_archive_joins': archive_inputs,
    'unchanged_original_apparatus': P['apparatus'], 'resources': len(resources), 'sequence': freeze['sequence'],
    'scope': 'Exact accepted producer/input reuse, no recompilation or historical performance samples. Original Create/all11 numerical API windows retain included candidate costs; fixture device negotiation remains outside those windows. Full owned wall retained.'})
print(json.dumps({'prepared': True, 'sequence': freeze['sequence'], 'resources': len(resources), 'executables': executables}))
