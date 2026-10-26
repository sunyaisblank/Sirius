"""Compose two verified Endpoint routing objects into the existing AO provider."""
from pathlib import Path
import difflib
import hashlib
import json
import subprocess
import sys
from component_gate import component, contract_inputs, doc, pin

R = Path.cwd()
W = Path(__file__).resolve().parent
P = W / 'provider'
B = R / 'bin/linux-gcc/native-cpu-endpoint-header-provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
CACHE = R / 'bin/linux-gcc/lavapipe-acqrel-ordering'
OLD = R / 'out/native-cpu-v4-provider/provider'
D = R / 'out/native-cpu-endpoint-header-domain'


def data_pin(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def canonical(value):
    return data_pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())


def write(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')


phase, = sys.argv[1:]
if phase == 'derive':
    replacements = {Path(row['path']).name: row for row in component()}
    assert not B.exists() and not P.exists()
    old = doc(OLD / 'source-contract.json')
    contract_inputs(old)
    assert doc(OLD / 'linked-native-proof.json')['passed']
    before = {str(path.relative_to(S)): pin(path) for path in S.rglob('*') if path.is_file()}
    assert canonical(before) == old['full_candidate_source_identity']
    objects = list(old['native_objects'])
    assert len(objects) == 8
    objects[1], objects[3] = replacements['endpoint.o'], replacements['admission.o']
    for row in objects + old['adapter_input_objects']:
        assert pin(R / row['path']) == {key: row[key] for key in ('bytes', 'sha256')}
    name = 'src/gallium/targets/lavapipe/meson.build'
    text = (S / name).read_text()
    for index in (1, 3):
        a, b = str(R / old['native_objects'][index]['path']), str(R / objects[index]['path'])
        assert text.count(a) == 1
        text = text.replace(a, b, 1)
    texts = {name: text}
    for name, a, b in (
        ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-v4-arithmetic-v39', 'Sirius-native-cpu-endpoint-header-v40'),
        ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8AO', 'SIR-CPU-25.2.8AP')):
        text = (S / name).read_text()
        assert text.count(a) == 1
        texts[name] = text.replace(a, b, 1)
    patch, after, modified = '', dict(before), {}
    for name, text in texts.items():
        patch += ''.join(difflib.unified_diff((S / name).read_text().splitlines(True),
            text.splitlines(True), n=3, fromfile='a/' + name, tofile='b/' + name))
        after[name] = data_pin(text.encode())
        modified[name] = dict(basis=before[name], candidate=after[name])
    inputs = dict(old['source_inputs'])
    for row in objects:
        inputs[row['path']] = {key: row[key] for key in ('bytes', 'sha256')}
    paths = [OLD / 'source-contract.json', OLD / 'linked-native-proof.json']
    paths += [path for path in D.iterdir() if path.is_file()]
    paths += [D / 'build/owner.json', D / 'compatibility/owner.json']
    paths += [path for path in W.iterdir() if path.is_file()]
    for path in paths:
        inputs[str(path.relative_to(R))] = pin(path)
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    P.mkdir()
    (P / 'provider.patch').write_text(patch)
    write(P / 'source-contract.json', dict(
        basis_source_identity=old['full_candidate_source_identity'], full_candidate_source_identity=canonical(after),
        modified_files=modified, patch=data_pin(patch.encode()), native_objects=objects,
        adapter_input_objects=old['adapter_input_objects'], component_manifest=old['component_manifest'],
        source_inputs=inputs, new_native_compilations=2, reused_raw_objects=3,
        guard_shared_byte_authority='out/native-cpu-endpoint-header-domain/native_packet.hpp',
        header_domain_component=dict(derivation='out/native-cpu-endpoint-header-domain/derivation.json',
            compatibility='out/native-cpu-endpoint-header-domain/compatibility.json',
            code_proof='out/native-cpu-endpoint-header-domain/code-proof.json',
            native_object_indices=[1, 3], relative_program_word=11241),
        scope='AP replaces only Endpoint scalar routing and pure admission with the verified shared predicate admitting one runtime Endpoint header word11241. The unchanged generated shader owns its acceptance/refusal and all output behavior. Other six AO objects, raw3 and merged original-v3/optional-v4/dispatcher math are byte-preserved; no merged relink or new raw compilation. Original bounds/alias/transactional, all other program bytes, Transport/DopriPhase, physical/count/precision/FP/policy/scheduler/host/SPIRV contracts remain. Finite component/source/ELF evidence does not transfer full38/C8/control/benefit/frame/scientific/release acceptance.'))
    print(json.dumps(dict(provider_patch_files=3, new_native_compilations=2, original_raw_objects_reused=3,
        reused_native_objects=6, new_merged_links=0)))
elif phase == 'configure':
    component()
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    contract = doc(P / 'source-contract.json')
    contract_inputs(contract)
    assert doc(P / 'review-admission.json')['passed']
    old = doc(OLD / 'actions.json')['build']
    inputs = {name for name in old['frozen_inputs'] if not (R / name).is_relative_to(CACHE)}
    inputs |= set(contract['source_inputs'])
    inputs |= {str(path.relative_to(R)) for path in S.rglob('*') if path.is_file()}
    inputs |= {str(path.relative_to(R)) for path in W.rglob('*') if path.is_file()}
    mutable_dirs = {'meson-private', 'meson-logs', 'meson-info', 'meson-uninstalled'}
    mutable_names = {'build.ninja', 'compile_commands.json', '.ninja_log', '.ninja_deps', '.ninja_lock'}
    for path in CACHE.rglob('*'):
        if path.is_file() and path.relative_to(CACHE).parts[0] not in mutable_dirs and path.name not in mutable_names:
            inputs.add(str(path.relative_to(R)))
    assert all((R / name).is_file() for name in inputs)
    action = dict(old, revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='configure_only', argv=['/usr/bin/python3', '-B', str(W / 'configure_provider.py')],
        frozen_inputs=sorted(inputs), limits=dict(seconds=60, rss_kib=4194304, reserve_kib=2097152),
        scope='Reconfigure the same Mesa cache/options/toolchain for two sealed replacement Endpoint routing objects; other six objects and original merged math unchanged. No compilation/provider workload or qualification claim.')
    write(P / 'actions.json', {'configure': action})
    for name in ('run_owned.py', 'owned_guard.py'):
        (P / name).write_bytes((W / name).read_bytes())
    print(json.dumps(dict(configure_inputs=len(inputs), new_merged_links=0)))
else:
    raise ValueError(phase)
