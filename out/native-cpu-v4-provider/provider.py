"""Compose the verified optional math into the unchanged private provider."""
from pathlib import Path
import difflib
import hashlib
import json
import subprocess
import sys
import zipfile

R = Path.cwd()
W = Path(__file__).resolve().parent
P = W / 'provider'
B = R / 'bin/linux-gcc/native-cpu-v4-provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
CACHE = R / 'bin/linux-gcc/lavapipe-acqrel-ordering'
OLD = R / 'out/native-cpu-local-identities/provider'
COMPONENT = R / 'attestations/qualification/fecafd04/native-cpu-v4-arithmetic-finite-compatibility/MANIFEST.json'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def data_pin(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def canonical(value):
    return data_pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())


def write(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')


def component():
    value = doc(COMPONENT)
    archive = R / value['archive']['path']
    assert pin(archive) == {k: value['archive'][k] for k in ('bytes', 'sha256')}
    assert value['finite_component_passed'] and not value['provider_adopted']
    for name, expected in value['products'].items():
        assert pin(R / name) == expected, name
    with zipfile.ZipFile(archive) as zipped:
        for name in ('code-proof.json', 'runtime-target-joins.json', 'dense-compatibility.json'):
            member = 'out/native-cpu-v4-arithmetic/' + name
            data = zipped.read(member)
            assert data_pin(data) == value['members'][member]
        assert json.loads(zipped.read('out/native-cpu-v4-arithmetic/code-proof.json'))['passed']
    return value


phase, = sys.argv[1:]
if phase == 'merge':
    value = component()
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    assert not B.exists() and not (W / 'actions.json').exists()
    B.mkdir()
    original = doc(R / 'out/native-cpu-value-fact-stores/prepared-cpp.json')['original_adapter']
    assert pin(R / original['path']) == {k: original[k] for k in ('bytes', 'sha256')}
    parts = [original] + [dict(path=name, **expected) for name, expected in value['products'].items()
                         if Path(name).name in ('v3.o', 'optional.o', 'dispatch.o')]
    assert len(parts) == 4
    destination = B / 'rows_adapter.o'
    argv = ['/usr/bin/ld', '-r', '-o', str(destination)] + [str(R / row['path']) for row in parts]
    inputs = {row['path'] for row in parts} | {
        str(COMPONENT.relative_to(R)), value['archive']['path'],
        'out/native-cpu-value-fact-stores/prepared-cpp.json', '/usr/bin/ld'}
    inputs |= {str(path.relative_to(R)) for path in W.iterdir() if path.is_file()}
    write(W / 'merge-inputs.json', dict(parts=parts, argv=argv,
          component_manifest=dict(path=str(COMPONENT.relative_to(R)), **pin(COMPONENT))))
    inputs.add(str((W / 'merge-inputs.json').relative_to(R)))
    write(W / 'actions.json', {'merge': dict(
        revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='build', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=argv, frozen_inputs=sorted(inputs),
        scope='One relocatable link of unchanged row adapter and verified v3/optional-v4/baseline-dispatch objects. No compilation, provider loading or numerical/performance claim.')})
elif phase == 'derive':
    owner = doc(W / 'merge/owner.json')
    assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    value = component()
    merged = B / 'rows_adapter.o'
    assert merged.is_file()
    P.mkdir(exist_ok=False)
    old = doc(OLD / 'source-contract.json')
    before = {str(path.relative_to(S)): pin(path) for path in S.rglob('*') if path.is_file()}
    assert canonical(before) == old['full_candidate_source_identity']
    parts = doc(W / 'merge-inputs.json')['parts']
    objects = old['native_objects'][:7] + [dict(path=str(merged.relative_to(R)), **pin(merged))]
    for row in objects:
        assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}
    texts = {}
    name = 'src/gallium/targets/lavapipe/meson.build'
    text = (S / name).read_text()
    a, b = str(R / old['native_objects'][7]['path']), str(merged)
    assert text.count(a) == 1
    texts[name] = text.replace(a, b, 1)
    for name, a, b in (
        ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-local-identities-v38', 'Sirius-native-cpu-v4-arithmetic-v39'),
        ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8AN', 'SIR-CPU-25.2.8AO')):
        text = (S / name).read_text()
        assert text.count(a) == 1
        texts[name] = text.replace(a, b, 1)
    patch, after, modified = '', dict(before), {}
    for name, text in texts.items():
        patch += ''.join(difflib.unified_diff((S / name).read_text().splitlines(True), text.splitlines(True), n=3, fromfile='a/' + name, tofile='b/' + name))
        after[name] = data_pin(text.encode())
        modified[name] = dict(basis=before[name], candidate=after[name])
    inputs = dict(old['source_inputs'])
    for row in parts + objects:
        inputs[row['path']] = {k: row[k] for k in ('bytes', 'sha256')}
    for name in value['products']:
        inputs[name] = value['products'][name]
    for path in (COMPONENT, R / value['archive']['path'], W / 'merge-inputs.json', W / 'merge/owner.json', Path(__file__)):
        inputs[str(path.relative_to(R))] = pin(path)
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    (P / 'provider.patch').write_text(patch)
    write(P / 'source-contract.json', dict(
        basis_source_identity=old['full_candidate_source_identity'], full_candidate_source_identity=canonical(after),
        modified_files=modified, patch=data_pin(patch.encode()), native_objects=objects,
        adapter_input_objects=parts, component_manifest=str(COMPONENT.relative_to(R)),
        source_inputs=inputs, new_native_compilations=0, reused_raw_objects=3,
        guard_shared_byte_authority=old['guard_shared_byte_authority'],
        scope='AO replaces only the merged arithmetic object with exact original-v3 fallback, optional same-eight-lane v4 math and baseline CPU dispatch. First7 native objects, raw3/scalar tails, packet/tuple/precision/physical/transactional contracts and original FP scopes/scheduler/host/SPIRVs remain. Provider identity strings distinguish this unadopted private candidate. New optional math fixture/code/actual-dense-PC evidence is source-bound in its pushed component bundle; no full38/C8/benefit/frame/scientific/release acceptance is transferred.'))
    print(json.dumps(dict(provider_patch_files=3, new_native_compilations=0, original_raw_objects_reused=3)))
elif phase == 'configure':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    assert doc(P / 'review-admission.json')['passed']
    old = doc(OLD / 'actions.json')['build']
    retired_path = R / 'attestations/qualification/86843614/native-cpu-local-identities-finite-compatibility/RETIRED-EXPANDED-LEDGERS.json'
    retired = doc(retired_path)
    assert retired['status'] == 'complete' and all(row['removed'] for row in retired['members'])
    excluded = {row['path'] for row in retired['members']}
    inputs = {name for name in old['frozen_inputs'] if not (R / name).is_relative_to(CACHE) and name not in excluded}
    inputs |= set(doc(P / 'source-contract.json')['source_inputs'])
    inputs |= {str(path.relative_to(R)) for path in S.rglob('*') if path.is_file()}
    inputs |= {str(path.relative_to(R)) for path in W.rglob('*') if path.is_file()}
    mutable_dirs = {'meson-private', 'meson-logs', 'meson-info', 'meson-uninstalled'}
    mutable_names = {'build.ninja', 'compile_commands.json', '.ninja_log', '.ninja_deps', '.ninja_lock'}
    for path in CACHE.rglob('*'):
        if path.is_file() and path.relative_to(CACHE).parts[0] not in mutable_dirs and path.name not in mutable_names:
            inputs.add(str(path.relative_to(R)))
    inputs.add(str(retired_path.relative_to(R)))
    assert all((R / name).is_file() for name in inputs)
    action = dict(old, revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='configure_only', argv=['/usr/bin/python3', '-B', str(W / 'configure_provider.py')],
        frozen_inputs=sorted(inputs), limits=dict(seconds=60, rss_kib=4194304, reserve_kib=2097152),
        scope='Reconfigure the existing unchanged Mesa cache/options/toolchain for one merged optional arithmetic object; no compiler/native workload or qualification claim.')
    write(P / 'actions.json', {'configure': action})
    for name in ('run_owned.py', 'owned_guard.py'):
        (P / name).write_bytes((W / name).read_bytes())
    print(json.dumps(dict(configure_inputs=len(inputs), excluded_verified_retired_ledgers=sorted(excluded.intersection(old['frozen_inputs'])))))
else:
    raise ValueError(phase)
