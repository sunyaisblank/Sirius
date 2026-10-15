"""Bounded original default-profile products, producer and eight-control owner actions."""
from pathlib import Path
import ast, hashlib, json, re, shlex, shutil, subprocess, sys
sys.dont_write_bytecode = True
W = Path(__file__).resolve().parent
R = W.parents[1]
T = R / 'bin/linux-gcc/portable-divide-fp32-estimate-restoring-controls'
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

def verify_host_reuse():
    from prerequisite_receipts import verify_canonical_producer
    return verify_canonical_producer(P,R)

def verify_output_seal(label):
    owner=accepted(label); join_runtime_inputs(label,owner)
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

def verify_prerequisites():
    from prerequisite_receipts import verify
    return verify(P,R)


def validate_actions():
    actions=json.loads((W/'actions.json').read_text())
    expected=[*P['producer_labels'],*[f'control-{i}' for i in range(1,9)],'verify-controls','verify-endpoint-routes']
    assert set(actions)==set(expected) and len(expected)==15
    required={'kind','revision','argv','frozen_inputs','limits','scope','expected_exit','software_provider'}
    outputs={name for names in P['outputs'].values() for name in names}
    for label in expected:
        a=actions[label];test=label.startswith('control-')
        assert set(a)==required and type(a['expected_exit']) is int and a['expected_exit']==0,label
        assert a['revision']==P['runner_preparation_revision'] and a['kind']==('test' if test else 'diagnostic')
        assert a['software_provider'] is test and isinstance(a['scope'],str) and a['scope']
        assert a['limits']=={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
        assert all(type(v) is int for v in a['limits'].values())
        assert a['argv']==(P['control_argv'][label] if test else ['/usr/bin/python3','-B',str(W/'gate.py'),'phase',label])
        files=a['frozen_inputs'];assert isinstance(files,list) and all(isinstance(n,str) and Path(n).is_absolute() for n in files)
        assert len(files)==len(set(files))
        assert not set(files)&set(P['outputs'].get(label,[]))
        assert not {str(W/label/n) for n in ['owner.json','inputs-before.json','inputs-after.json','output-seal.json']}&set(files)
        assert not {str(W/n) for n in P['outside_inputs']}&set(files)
    assert (W/'actions.json').read_bytes()==(W/'actions.unpinned.json').read_bytes()
    return actions

def pin_context():
    validate_actions(); verify_static(); verify_host_reuse(); verify_prerequisites()
    assert not T.exists() and not (W / 'operator-context.json').exists()
    assert not subprocess.check_output(['git', 'status', '--porcelain=v1'], cwd=R)
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=R, text=True).strip()
    assert (W / 'actions.json').read_bytes() == (W / 'actions.unpinned.json').read_bytes()
    assert revision == P['runner_preparation_revision']
    actions = json.loads((W / 'actions.json').read_text())
    assert all(action['revision'] == revision for action in actions.values())
    write(W / 'operator-context.json', {'revision': revision, 'actions': record(W / 'actions.json'),
          'source_freeze': record(W / 'source-freeze.json'), 'preparation': record(W / 'preparation.json')})
    print('Pinned restoring-FP32 controls source; no compiler/provider launched.', flush=True)

def prepare_producer():
    verify_prerequisites(); verify_host_reuse()
    directory = T / 'controls'; assert not directory.exists()
    (directory / 'headers').mkdir(parents=True)
    for name in ('retained_kernels.h', 'retained_owner_kernels.h', 'retained_owner_layout.h'):
        shutil.copyfile(Path(P['candidate_products']) / name, directory / 'headers' / name)
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
            'scope': 'Actual candidate linked to exact accepted canonical Vulkan object/base. No private Int64/Float64 device TU or paired55-feature requirement and no fresh numerical feature-argument trace is claimed. Original baseline ELF remains unchanged; no production adoption or qualification.'}

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
    validate_actions(); verify_static()
    context = json.loads((W / 'operator-context.json').read_text())
    assert context['actions'] == record(W / 'actions.json') and context['preparation'] == record(W / 'preparation.json')
    assert context['source_freeze'] == record(W / 'source-freeze.json')
    for prior in ancestor_labels(label): verify_output_seal(prior)
    assert all(not Path(path).exists() for path in P['outputs'][label])
    if label == 'prepare-producer': prepare_producer()
    else:
        for argv in P['recipes'][label]: subprocess.run(argv, cwd=R, check=True)
        if label == 'verify-producer': write(W / 'producer-feature-acceptance.json', producer_feature_review())
    write(W / label / 'output-seal.json', {'label': label, 'recipe': P['recipes'][label],
          'source_freeze': record(W / 'source-freeze.json'), 'operator_context': record(W / 'operator-context.json'),
          'inputs_before': record(W / label / 'inputs-before.json'),
          'outputs': {name: record(name) for name in P['outputs'][label]}})

def pin_runtime():
    validate_actions(); verify_static(); prerequisites=verify_prerequisites(); verify_host_reuse()
    assert not (W / 'runtime-context.json').exists()
    for label in P['producer_labels']: verify_output_seal(label)
    assert json.loads((W / 'producer-acceptance.json').read_text())['passed']
    assert json.loads((W / 'producer-feature-acceptance.json').read_text())['passed']
    write(W / 'runtime-context.json', {'profile': record(W / 'profile.json'),
          'producer': record(W / 'producer-acceptance.json'), 'feature': record(W / 'producer-feature-acceptance.json'),
          'inherited_full_profile': prerequisites, 'arithmetic': prerequisites['arithmetic'],
          'producer_output_seals': {label: verify_output_seal(label) for label in P['producer_labels']},
          'executable': record(T / 'controls/original_candidate_controls'),
          'scope': 'Original eight controls only, following5fresh producer owners and immutable2producer+1separateprofile receipts. All original filters/assertions/limits/provider/cache policy retained.'})
    print('Original eight controls admitted behind terminal compiler/profile/producer acceptance.', flush=True)

def review():
    validate_actions(); verify_static(); prerequisites=verify_prerequisites(); verify_host_reuse()
    runtime = json.loads((W / 'runtime-context.json').read_text())
    assert runtime['inherited_full_profile'] == prerequisites
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
    write(W / 'prerequisites-result.json', {'schema':'sirius-fp32-restoring-original-controls-result-v1','passed': True, 'original_eight_controls_passed': True,
          'profile': record(W / 'profile.json'), 'arithmetic': prerequisites['arithmetic'],
          'control_acceptance': record(W / 'control-acceptance.json'), 'endpoint_route_acceptance': record(W / 'endpoint-route-acceptance.json'),
          'production_adopted': False, 'benefit_executed': False,
          'fresh_owner_count':15,'fresh_producer_count':5,'inherited_prerequisites':prerequisites,
          'scope': 'Original primitive arithmetic and eight integrated prerequisites for explicit17-module guarded FP32 restoring diagnostic profile. Original20+2 adoption gates remain unexecuted; no full-science/frame/release qualification.'})
    print('Original eight controls and complete route/readback joins accepted; no benefit verdict.', flush=True)

if __name__ == '__main__':
    if sys.argv[1] == 'phase': phase(sys.argv[2])
    else: {'pin-context': pin_context, 'pin-runtime': pin_runtime, 'review': review}[sys.argv[1]]()
