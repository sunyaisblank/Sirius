# Derived from the exercised S frame planner. Only current provider/component identities change; original resources/profile/governors retained.
"""Original frame with fresh AB full38/C8 and explicitly scoped Z predecessor gates."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
P = W.parent
Q = R / 'out/native-cpu-direct-scratch-cohort'
F = R / 'out/native-cpu-dispatch-admission/cold-frame'
MANIFEST = R / 'attestations/qualification/629e3794/native-cpu-dispatch-admission-frame-timeout/MANIFEST.json'
EXE = R / 'bin/linux-gcc/native-cpu-upload-aware-staging/frame/sirius_render_tests'
ARCHIVE = R / 'bin/linux-gcc/native-cpu-upload-aware-staging/matched/libsirius_backend.a'
HOST = ARCHIVE.parent / 'audited_backend_tests'
LIMITS = {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
controls = doc(Q / 'controls/readers/controls-1-8-acceptance.json')
assert controls['all8_original_controls_pass'] and controls['all22_original_benefit_checks_pass']
benefit = doc(Q / 'benefit/readers/benefit.json')
decimal = doc(Q / 'benefit/readers/decimal-agreement.json')
assert benefit['audit_passed'] and benefit['all22_numeric_pass'] and benefit['checks_passed'] == 22
assert decimal['agrees_on_all22'] and decimal['all22_pass']
c8 = doc(P / 'control-8-acceptance.json')
assert c8['passed'] and c8['all32_readbacks_match_preserved_originals']
assert c8['native_first_packets_handled'] == c8['pure_preflights_handled'] == 3
expected_raw = (2 + int(any(row['selected'] and row['rows'] >= 8 for row in c8['observed_transport_policy'])) if c8['runtime_SIMD_supported'] else 0)
assert c8['first_raw_SIMD_entries_completed'] == expected_raw
assert c8['observed_Transport12_implementation'] in ('F scalar one-row tasks', 'SIMD8 plus scalar tails')
assert c8['transport_component_scope']['current_AB_link_proof'] == pin(R / 'out/native-cpu-predicated-products/provider/linked-native-proof.json')
assert doc(R / 'out/native-cpu-predicated-products/provider/linked-native-proof.json')['passed']
provider = doc(R / 'out/native-cpu-predicated-products/provider/provider-build-result.json')
assert provider['passed']
provider = {**provider, **{key: str(R / provider[key]) for key in ('library', 'ICD')}}
assert not doc(Q / 'cold-frame/result.json')['accepted']
assert doc(R / 'out/native-cpu-predicated-products/ispc-derivation.json')['instructions'] == 9380
for key in ('library', 'ICD'):
    assert pin(Path(provider[key])) == provider[key + '_identity']

assert pin(MANIFEST)['sha256'] == 'e82c6e78132fbb34533a8d284f949bda60f78e0d23773861ef3b6478eadcbd69'
assert pin(F / 'build-result.json')['sha256'] == '558c98ca9463708223e20ba0987f8426f803165757a52c266350d736da4284f3'
assert pin(F / 'link-recipe.json')['sha256'] == '8cba87eb7c419379069e22623dace65d80601abaf1497598c3bf8b232c790e31'
assert pin(F / 'actions.json')['sha256'] == 'ed6f4d3daa35f228ad15887a7426f9c2c3f03befd8ac79d642715688fb46d733'
build = doc(F / 'build-result.json')
recipe = doc(F / 'link-recipe.json')
assert build['passed'] and build['original_objects_reused'] and build['new_CPP_compilations'] == 0
host = doc(R / 'out/native-cpu-upload-aware-staging/host-build-result.json')
assert host['passed'] and host['new_cpp_compilations'] == 1
assert host['all3modules_once_in_both_executables']
assert str(EXE) == str(R / host['frame_executable'])
assert pin(EXE) == host['frame_executable_identity']
assert str(ARCHIVE) == str(R / host['backend_archive'])
assert pin(ARCHIVE) == host['backend_archive_identity']
assert str(HOST) == str(R / host['executable'])
assert pin(HOST) == host['executable_identity']
assert pin(R / 'tests/operational/vulkan_runtime_profile.cmake') == {
    'bytes': 990, 'sha256': 'e39e8aee8f82b4f1934a42542cd5b743822bed02f6ed15e1a3eaa6185bc0f8ea'}
assert recipe['argv'].count(str(R / 'bin/linux-gcc/native-cpu-row-ownership/matched/libsirius_backend.a')) == 1
assert recipe['argv'][recipe['argv'].index('-o') + 1] == str(R / 'bin/linux-gcc/native-cpu-dispatch-admission/frame/sirius_render_tests')
recovery = doc(MANIFEST)['resources_recovery']
cache, destination = R / recovery['source'], EXE.parent / 'resources'
assert destination == EXE.parent / 'resources'
expected = {name: {k: value[k] for k in ('bytes', 'sha256')}
            for name, value in recovery['files'].items()}
assert len(expected) == 32 and sum(value['bytes'] for value in expected.values()) == 54156187
assert expected == recipe['resources'] == build['original_resources_copied_exactly']
for name, value in expected.items():
    relative = Path(name)
    assert not relative.is_absolute() and '..' not in relative.parts
    source = cache / relative
    assert source.is_file() and not source.is_symlink() and source.resolve().is_relative_to(cache.resolve())
    assert pin(source) == value
    assert source.stat().st_mode & 0o777 == recovery['files'][name]['mode']
    assert source.stat().st_mtime_ns == recovery['files'][name]['mtime_ns']
modules = doc(R / 'out/native-cpu-row-ownership/compiler-result.json')['modules']
binary = EXE.read_bytes()
joins = {}
for kind in ('Transport', 'Endpoint', 'DopriPhase'):
    module = modules[kind]
    path = R / module['module']
    assert pin(path) == {k: module[k] for k in ('bytes', 'sha256')}
    assert binary.count(path.read_bytes()) == 1
    joins[kind] = {'module': module['module'], **pin(path), 'complete_occurrences': 1}
assert joins == build['module_joins']
del binary

mode, = sys.argv[1:]
if mode == 'restore':
    assert not destination.exists() and not destination.is_symlink()
    assert not (W / 'resource-restore.json').exists() and not (W / 'reused-frame.json').exists()
    assert cache.stat().st_dev == destination.parent.stat().st_dev
    journal = {'status': 'started', 'source': revision, 'manifest': pin(MANIFEST),
               'cache': str(cache.relative_to(R)), 'destination': str(destination.relative_to(R)),
               'method': 'hardlinks to frozen original cache; no new resource payload', 'files': {}}
    with (W / 'resource-restore.json').open('x') as stream:
        json.dump(journal, stream, indent=2)
        stream.write('\n')
    destination.mkdir()
    for name, value in expected.items():
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        os.link(cache / name, target)
        assert os.path.samefile(cache / name, target) and pin(target) == value
        journal['files'][name] = value
        write(W / 'resource-restore.json', journal)
    assert {str(path.relative_to(destination)) for path in destination.rglob('*') if path.is_file()} == set(expected)
    journal.update(status='complete', resource_payload_bytes_copied=0)
    write(W / 'resource-restore.json', journal)
    write(W / 'frame-inputs.json', {'passed': True, 'source': revision,
        'executable': str(EXE.relative_to(R)), 'executable_identity': pin(EXE),
        'historical_build': pin(F / 'build-result.json'), 'historical_link_recipe': pin(F / 'link-recipe.json'),
        'host_build': pin(R / 'out/native-cpu-upload-aware-staging/host-build-result.json'),
        'backend_archive': pin(ARCHIVE), 'matched_host': pin(HOST), 'module_joins': joins,
        'resource_restore': pin(W / 'resource-restore.json'), 'new_links': 2, 'new_CPP_compilations': 1})
    print(json.dumps({'restored_names': 32, 'resource_payload_bytes_copied': 0, 'host_links_already_completed': 2, 'executed': False}))
elif mode == 'frame-actions':
    assert not (W / 'actions.json').exists() and not (W / 'frame').exists()
    restored = doc(W / 'resource-restore.json')
    reuse = doc(W / 'frame-inputs.json')
    assert restored['status'] == 'complete' and restored['files'] == expected
    assert reuse['passed'] and reuse['new_links'] == 2 and reuse['new_CPP_compilations'] == 1
    assert reuse['host_build'] == pin(R / 'out/native-cpu-upload-aware-staging/host-build-result.json')
    assert reuse['executable_identity'] == pin(EXE) and reuse['module_joins'] == joins
    assert reuse['resource_restore'] == pin(W / 'resource-restore.json')
    assert {str(path.relative_to(destination)) for path in destination.rglob('*') if path.is_file()} == set(expected)
    for name, value in expected.items():
        assert os.path.samefile(cache / name, destination / name) and pin(destination / name) == value
    dependencies = set()
    for path in (EXE, Path('/usr/bin/cmake'), Path(provider['library'])):
        text = subprocess.check_output(['ldd', str(path)], text=True)
        dependencies |= {Path(value).resolve() for value in re.findall(r'(/[^\s()]+)', text) if Path(value).is_file()}
    c8_inputs = doc(P / 'control-8/inputs-before.json')
    for path in dependencies:
        if str(path) in c8_inputs:
            assert pin(path) == {key: c8_inputs[str(path)][key] for key in ('bytes', 'sha256')}
    inputs = {R / name for name in doc(Q / 'controls/runtime/actions.json')['control-7']['frozen_inputs']}
    inputs |= {R / name for name in doc(P / 'actions.json')['control-8']['frozen_inputs']}
    inputs |= dependencies | {EXE, ARCHIVE, HOST, MANIFEST, F / 'build-result.json', F / 'link-recipe.json', F / 'actions.json',
        Path(provider['library']), Path(provider['ICD']), Path('/usr/bin/cmake'), Path('/usr/bin/python3'), Path('/etc/ld.so.cache'),
        R / 'tests/operational/vulkan_runtime_profile.cmake', R / 'out/native-cpu-row-ownership/compiler-result.json',
        R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py', R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py',
        R / 'out/native-cpu-dispatch-admission/cold-frame/read_result.py', R / 'out/native-cpu-arithmetic-sharing/readers/native_evidence.py'}
    inputs |= {path for path in W.iterdir() if path.is_file()}
    inputs |= {R / 'out/native-cpu-upload-aware-staging/host-build-result.json',
               R / 'out/native-cpu-predicated-products/provider/linked-native-proof.json',
               R / 'out/native-cpu-runtime-observe/all_threads.py'}
    inputs |= {cache / name for name in expected} | {destination / name for name in expected}
    inputs |= {R / 'out/native-cpu-predicated-products/provider/provider-build-result.json',
               R / 'out/native-cpu-predicated-products/ispc-derivation.json',
               Q / 'cold-frame/result.json'}
    for name in ('controls/readers/controls-1-8-acceptance.json', 'benefit/readers/benefit.json',
                 'benefit/readers/decimal-agreement.json', 'control-8-acceptance.json',
                 'control-8/native-observation.json', 'control-8/inputs-before.json'):
        inputs.add((Q if name.startswith(('controls/', 'benefit/')) else P) / name)
    for control in range(1, 8):
        for name in ('owner.json', 'gtest.xml', 'inputs-before.json', 'inputs-after.json'):
            inputs.add(Q / 'controls/runtime' / ('control-' + str(control)) / name)
    assert all(path.is_file() for path in inputs)
    temp = R / ('renders/native-cpu-predicated-products-cohort-' + revision[:8] + '/cold-frame-tmp')
    assert not temp.exists()
    temp.mkdir(parents=True)
    action = doc(F / 'actions.json')['frame']
    assert action['limits'] == LIMITS
    assert action['argv'] == ['/usr/bin/cmake', '-DSIRIUS_RENDER_TESTS=' + str(R / 'bin/linux-gcc/native-cpu-dispatch-admission/frame/sirius_render_tests'), '-P',
                              str(R / 'tests/operational/vulkan_runtime_profile.cmake')]
    action['argv'][1] = '-DSIRIUS_RENDER_TESTS=' + str(EXE)
    action.update(revision=revision, frozen_inputs=sorted(str(path.relative_to(R)) if path.is_relative_to(R) else str(path) for path in inputs),
                  runtime_provider=provider['library'], runtime_icd=provider['ICD'], build_environment={'TMPDIR': str(temp)},
                  scope=action['scope'] + ' Genuine new host/frame ELF with one changed C++ unit and two audited links; exact AB provider with predicated product accumulators; original complete nonzero cores and invalid/zero results retained, AA metadata/scratch ABI and Y scheduler/FP scope retained. Fresh AB full38/C8; Z all22+Decimal/all8 are predecessor evidence only, not newly executed AA gates. Only the32 byte-verified historical resource names restored as hardlinks.')
    write(W / 'actions.json', {'frame': action})
    print(json.dumps({'prepared_original_frame': True, 'host_links_already_completed': 2, 'executed': False}))
elif mode == 'retire-resources':
    assert not (W / 'resource-retirement.json').exists()
    result = doc(W / 'result.json')
    owner = doc(W / 'frame/owner.json')
    assert result['attempted'] and result['owner'] == pin(W / 'frame/owner.json')
    assert owner['status'] == 'terminal' and owner['remaining_owned_processes'] == []
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and not owner['cleanup_errors']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    restored = doc(W / 'resource-restore.json')
    assert restored['status'] == 'complete' and restored['files'] == expected
    assert {str(path.relative_to(destination)) for path in destination.rglob('*') if path.is_file()} == set(expected)
    for name, value in expected.items():
        assert os.path.samefile(cache / name, destination / name)
        assert pin(cache / name) == pin(destination / name) == value
    retirement = {'status': 'started', 'result': pin(W / 'result.json'),
                  'verified_original_cache_and_restored_files': expected, 'unlinked_restored_names': []}
    write(W / 'resource-retirement.json', retirement)
    directories = {destination}
    for name in expected:
        target = destination / name
        parent = target.parent
        while parent != destination:
            directories.add(parent)
            parent = parent.parent
        target.unlink()  # Remove only this verified additional name.
        retirement['unlinked_restored_names'].append(name)
        write(W / 'resource-retirement.json', retirement)
    for directory in sorted(directories, key=lambda path: len(path.parts), reverse=True):
        directory.rmdir()  # Fail on any unrecorded file; never remove a tree.
    assert all(pin(cache / name) == value for name, value in expected.items())
    retirement.update(status='complete', original_cache_preserved=True)
    write(W / 'resource-retirement.json', retirement)
    print(json.dumps({'unlinked_restored_names': 32, 'original_cache_preserved': True}))
else:
    raise ValueError(mode)
