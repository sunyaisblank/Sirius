"""Prepare AV's three-file patch using the unchanged AU headers/cache."""
from pathlib import Path
import difflib
import hashlib
import json
import subprocess
import sys
from component_gate import component, contract_inputs, doc, pin, HEADERS

R = Path.cwd()
W = Path(__file__).resolve().parent
P = W / 'provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
CACHE = R / 'bin/linux-gcc/lavapipe-acqrel-ordering'
OLD = R / 'out/native-cpu-shared-value-provider/provider'
IDENTITY = 'Sirius-native-cpu-closed-dp-v46'
UUID = 'SIR-CPU-25.2.8AV'


def data_pin(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def canonical(value):
    return data_pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


phase, = sys.argv[1:]
value = component()
old = doc(OLD / 'source-contract.json')
if phase == 'derive':
    assert not P.exists()
    assert doc(OLD / 'linked-native-proof.json')['passed']
    # Identity collision is checked against existing provider patches and Mesa.
    for path in (R / 'out').glob('native-cpu-*-provider/provider/provider.patch'):
        assert IDENTITY not in path.read_text() and UUID not in path.read_text(), path
    before = {str(path.relative_to(S)): pin(path) for path in S.rglob('*') if path.is_file()}
    assert canonical(before) == old['full_candidate_source_identity']
    objects = value['native_first7'] + [value['merged']]
    name = 'src/gallium/targets/lavapipe/meson.build'
    text = (S / name).read_text()
    a, b = str(R / old['native_objects'][7]['path']), str(R / objects[7]['path'])
    assert text.count(a) == 1 and b not in text
    texts = {name: text.replace(a, b, 1)}
    for name, a, b in (
        ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-shared-value-v45', IDENTITY),
        ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8AU', UUID)):
        text = (S / name).read_text()
        assert text.count(a) == 1 and b not in text
        texts[name] = text.replace(a, b, 1)
    patch, after, modified = '', dict(before), {}
    for name, text in texts.items():
        patch += ''.join(difflib.unified_diff((S / name).read_text().splitlines(True), text.splitlines(True),
            n=3, fromfile='a/' + name, tofile='b/' + name))
        after[name] = data_pin(text.encode())
        modified[name] = dict(basis=before[name], candidate=after[name])
    assert len(modified) == 3 and not set(modified).intersection(HEADERS)
    inputs = dict(old['source_inputs'])
    inputs.update(value['source_inputs'])
    paths = [OLD / 'source-contract.json', OLD / 'linked-native-proof.json', OLD / 'incremental-plan.json']
    paths += [path for path in W.iterdir() if path.is_file()]
    for path in paths:
        inputs[str(path.relative_to(R))] = pin(path)
    for row in objects + value['parts']:
        inputs[row['path']] = {key: row[key] for key in ('bytes', 'sha256')}
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    P.mkdir()
    (P / 'provider.patch').write_text(patch)
    write(P / 'source-contract.json', dict(
        basis_source_identity=old['full_candidate_source_identity'], full_candidate_source_identity=canonical(after),
        modified_files=modified, patch=data_pin(patch.encode()), native_objects=objects,
        adapter_input_objects=value['parts'], raw_v4_objects=value['raw_v4_objects'],
        shared_value_math_object=value['shared_value_math_object'], shared_value_interface=value['shared_value_interface'],
        closed_dp_math_object=value['closed_dp_math_object'], closed_dp_interface=value['closed_dp_interface'], source_inputs=inputs,
        component_interface=value['component_interface'], component_build_receipt='out/native-cpu-closed-dp-component/build-result.json',
        unchanged_bridge_headers=value['unchanged_bridge_headers'], predecessor_source_contract=str((OLD / 'source-contract.json').relative_to(R)),
        component_manifest=old['component_manifest'], admission_component_manifest=old['admission_component_manifest'],
        header_domain_component=old['header_domain_component'],
        guard_shared_byte_authority=old['guard_shared_byte_authority'],
        admission_header_patch_identity=old['admission_header_patch_identity'],
        observer_layout_assertion_identity=old['observer_layout_assertion_identity'],
        new_native_compilations=0, component_baseline_compilations=1, reused_native_indices=list(range(7)),
        reused_baseline_raw_objects=3, reused_optional_raw_objects=3,
        shared_value_raw_objects=2, shared_value_core_objects=1, closed_dp_raw_objects=1, closed_dp_core_objects=1,
        provider_identity=IDENTITY, provider_UUID=UUID,
        scope='AV changes only the eighth object path and two identities. First7 objects, headers, cs_exec_fn, all public/scalar/ISA/fallback branches and AU admitted T/E remain exact. Only admitted DP v4 selects the closed DP and existing local-effects strict cc92 math9. Twelve constituents include six old public AS/GCC/dispatcher objects, adapter, AU T/E/shared9, closed DP/local9. No scene/count/layout/allocation/scheduler/cache/governor/FP changes or numerical transfer. Fresh static source/PC/ISA/FP proof and original genuine-job C8/all32 are required. No ABBA/frame until remaining-stage closure is justified.'))
    print(json.dumps(dict(provider_patch_files=3, unchanged_bridge_headers=2, native_objects=8,
        merged_parts=12, numerical_qualification_transferred=False)))
elif phase == 'configure':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    contract = doc(P / 'source-contract.json')
    contract_inputs(contract)
    assert doc(P / 'review-admission.json')['passed']
    old_action = doc(OLD / 'actions.json')['build']
    inputs = {name for name in old_action['frozen_inputs'] if not (R / name).is_relative_to(CACHE)}
    inputs |= set(contract['source_inputs'])
    inputs |= {str(path.relative_to(R)) for path in S.rglob('*') if path.is_file()}
    inputs |= {str(path.relative_to(R)) for path in W.rglob('*') if path.is_file()}
    mutable_dirs = {'meson-private', 'meson-logs', 'meson-info', 'meson-uninstalled'}
    mutable_names = {'build.ninja', 'compile_commands.json', '.ninja_log', '.ninja_deps', '.ninja_lock'}
    for path in CACHE.rglob('*'):
        if path.is_file() and path.relative_to(CACHE).parts[0] not in mutable_dirs and path.name not in mutable_names:
            inputs.add(str(path.relative_to(R)))
    assert all((R / name).is_file() for name in inputs)
    action = dict(old_action, revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='configure_only', argv=['/usr/bin/python3', '-B', str(W / 'configure_provider.py')],
        frozen_inputs=sorted(inputs), limits=dict(seconds=60, rss_kib=4194304, reserve_kib=2097152),
        scope='Same Mesa cache/options/toolchain, eight native rows with the new twelve-input eighth object; both bridge headers unchanged. Configuration only; no numerical/FP/ISA/performance acceptance.')
    write(P / 'actions.json', {'configure': action})
    for name in ('run_owned.py', 'owned_guard.py'):
        (P / name).write_bytes((W / name).read_bytes())
    print(json.dumps(dict(configure_inputs=len(inputs), merged_parts=12)))
else:
    raise ValueError(phase)
