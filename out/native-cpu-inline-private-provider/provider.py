"""Prepare one seven-input merge and the AS three-file private-provider patch."""
from pathlib import Path
import difflib
import hashlib
import json
import subprocess
import sys
from component_gate import component, contract_inputs, doc, pin, terminal

R = Path.cwd()
W = Path(__file__).resolve().parent
P = W / 'provider'
B = R / 'bin/linux-gcc/native-cpu-inline-private-provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
CACHE = R / 'bin/linux-gcc/lavapipe-acqrel-ordering'
OLD = R / 'out/native-cpu-packed-math-provider/provider'


def data_pin(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def canonical(value):
    return data_pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())


def write(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')


def component_inputs(value):
    paths = [value['manifest'], value['archive'], value['math_manifest'], value['math_archive']]
    paths += [R / row['path'] for row in value['parts']]
    return {str(path.relative_to(R)): pin(path) for path in paths}


phase, = sys.argv[1:]
if phase == 'merge':
    value = component()
    old = doc(OLD / 'source-contract.json')
    contract_inputs(old)
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    assert not B.exists() and not (W / 'actions.json').exists()
    parts = value['parts']
    assert len(parts) == 7
    destination = B / 'rows_adapter.o'
    argv = ['/usr/bin/ld', '-r', '-o', str(destination)] + [str(R / row['path']) for row in parts]
    inputs = set(old['source_inputs']) | set(component_inputs(value)) | {'/usr/bin/ld'}
    inputs |= {str(path.relative_to(R)) for path in W.iterdir() if path.is_file()}
    B.mkdir()
    write(W / 'merge-inputs.json', dict(parts=parts, argv=argv,
        raw_component_manifest=dict(path=str(value['manifest'].relative_to(R)), **pin(value['manifest'])),
        raw_component_archive=dict(path=str(value['archive'].relative_to(R)), **pin(value['archive']))))
    inputs.add(str((W / 'merge-inputs.json').relative_to(R)))
    assert all((R / name).is_file() for name in inputs)
    write(W / 'actions.json', {'merge': dict(
        revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='build', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152), argv=argv,
        frozen_inputs=sorted(inputs),
        scope='One relocatable seven-input link: unchanged AQ baseline selector, sealed opcode-inline/privateD-zero raw3, unchanged AO baseline dispatcher/original-v3 and separately qualified zero-residual v4 math. No compilation, provider loading or numerical/performance claim.')})
elif phase == 'derive':
    terminal(W / 'merge/owner.json')
    value = component()
    old = doc(OLD / 'source-contract.json')
    contract_inputs(old)
    assert doc(OLD / 'linked-native-proof.json')['passed']
    merged = B / 'rows_adapter.o'
    assert merged.is_file() and not P.exists()
    parts = doc(W / 'merge-inputs.json')['parts']
    assert parts == value['parts'] and len(parts) == 7
    before = {str(path.relative_to(S)): pin(path) for path in S.rglob('*') if path.is_file()}
    assert canonical(before) == old['full_candidate_source_identity']
    objects = old['native_objects'][:7] + [dict(path=str(merged.relative_to(R)), **pin(merged))]
    for row in objects + parts:
        assert pin(R / row['path']) == {key: row[key] for key in ('bytes', 'sha256')}
    name = 'src/gallium/targets/lavapipe/meson.build'
    text = (S / name).read_text()
    a, b = str(R / old['native_objects'][7]['path']), str(merged)
    assert text.count(a) == 1
    texts = {name: text.replace(a, b, 1)}
    for name, a, b in (
        ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-packed-math-v42', 'Sirius-native-cpu-inline-private-v43'),
        ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8AR', 'SIR-CPU-25.2.8AS')):
        text = (S / name).read_text()
        assert text.count(a) == 1
        texts[name] = text.replace(a, b, 1)
    patch, after, modified = '', dict(before), {}
    for name, text in texts.items():
        patch += ''.join(difflib.unified_diff((S / name).read_text().splitlines(True), text.splitlines(True),
            n=3, fromfile='a/' + name, tofile='b/' + name))
        after[name] = data_pin(text.encode())
        modified[name] = dict(basis=before[name], candidate=after[name])
    inputs = dict(old['source_inputs'])
    inputs.update(component_inputs(value))
    for row in objects:
        inputs[row['path']] = {key: row[key] for key in ('bytes', 'sha256')}
    paths = [OLD / 'source-contract.json', OLD / 'linked-native-proof.json',
             W / 'merge-inputs.json', W / 'merge/owner.json']
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
        adapter_input_objects=parts, component_manifest=old['component_manifest'],
        raw_component_manifest=str(value['manifest'].relative_to(R)),
        raw_component_archive=str(value['archive'].relative_to(R)),
        math_component_manifest=str(value['math_manifest'].relative_to(R)),
        math_component_archive=str(value['math_archive'].relative_to(R)),
        source_inputs=inputs, new_native_compilations=0, sealed_component_native_compilations=3,
        reused_raw_objects=3, new_optional_raw_objects=3,
        guard_shared_byte_authority=old['guard_shared_byte_authority'],
        header_domain_component=old['header_domain_component'],
        scope='AS keeps first7 AP native objects and the original AQ baseline selector exact. The eighth replaces only AR raw3 with forced opcode-inline stages and private destination zeroing; the same baseline adapter, dispatcher/original-v3 and sealed zero-residual v4 math are reused. All tagged math calls use private576B temporary banks; remaining dense calls already use private local banks. All9380 ordered tuples remain; only the opcode function is forced inline, while generic fallback and399 layer definitions remain unchanged. Original fallback/admission, complete fields/RHS/trials/masks/FP/ISA contracts, scheduler/host/SPIRVs/counts/cache/governors remain. Static/finite component evidence does not transfer fresh C8/all8/benefit/frame/scientific/release acceptance; stack reuse, SSA and performance are unproven.'))
    print(json.dumps(dict(provider_patch_files=3, new_native_compilations=0, merged_parts=7,
        first7_native_objects_unchanged=True, sealed_component_native_compilations=3)))
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
        scope='Reconfigure the same Mesa cache/options/toolchain for one merged seven-input object; first7 AP objects unchanged. No compiler/native workload or qualification claim.')
    write(P / 'actions.json', {'configure': action})
    for name in ('run_owned.py', 'owned_guard.py'):
        (P / name).write_bytes((W / name).read_bytes())
    print(json.dumps(dict(configure_inputs=len(inputs), merged_parts=7)))
else:
    raise ValueError(phase)
