"""Prepare one seven-input merge and the AT five-file private-provider patch."""
from pathlib import Path
import difflib
import hashlib
import json
import subprocess
import sys
from component_gate import component, contract_inputs, doc, pin, terminal, header_patch, HEADERS, add_layout_assertions, remove_layout_assertions, LAYOUT_ASSERTIONS

R = Path.cwd()
W = Path(__file__).resolve().parent
P = W / 'provider'
B = R / 'bin/linux-gcc/native-cpu-dispatch-token-provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
CACHE = R / 'bin/linux-gcc/lavapipe-acqrel-ordering'
OLD = R / 'out/native-cpu-inline-private-provider/provider'


def data_pin(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def canonical(value):
    return data_pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())


def write(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')


def component_inputs(value):
    paths = [value['manifest'], value['archive'], value['math_manifest'], value['math_archive'],
             OLD / 'incremental-plan.json']
    paths += [R / row['path'] for row in value['parts'] + value['native_first7']]
    paths += [R / 'out/native-cpu-dispatch-admission-token' / name for name in value['source_derivation']['derived']]
    return {str(path.relative_to(R)): pin(path) for path in paths}


def invariant_basis_inputs(old, value):
    # Exactly these two declared patch targets cease being invariant inputs.
    # Their old identities remain explicit basis/inverse authorities.
    excluded={str(S/name):name for name in HEADERS}
    inputs={}
    for name,expected in old['source_inputs'].items():
        absolute=str(R/name)
        if absolute in excluded:
            assert expected==value['header_basis'][excluded[absolute]]
        else:
            assert pin(R/name)==expected,name
            inputs[name]=expected
    return inputs


phase, = sys.argv[1:]
if phase == 'merge':
    value = component()
    old = doc(OLD / 'source-contract.json')
    invariant_basis_inputs(old, value)
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
        scope='One relocatable seven-input link: new baseline public/private selector, unchanged sealed AS raw3, unchanged AO baseline dispatcher/original-v3 and separately qualified zero-residual v4 math. No compilation, provider loading or numerical/performance claim.')})
elif phase == 'derive':
    terminal(W / 'merge/owner.json')
    value = component()
    old = doc(OLD / 'source-contract.json')
    invariant_basis_inputs(old, value)
    assert doc(OLD / 'linked-native-proof.json')['passed']
    merged = B / 'rows_adapter.o'
    assert merged.is_file() and not P.exists()
    parts = doc(W / 'merge-inputs.json')['parts']
    assert parts == value['parts'] and len(parts) == 7
    before = {str(path.relative_to(S)): pin(path) for path in S.rglob('*') if path.is_file()}
    assert canonical(before) == old['full_candidate_source_identity']
    objects = value['native_first7'] + [dict(path=str(merged.relative_to(R)), **pin(merged))]
    for row in objects + parts:
        assert pin(R / row['path']) == {key: row[key] for key in ('bytes', 'sha256')}
    name = 'src/gallium/targets/lavapipe/meson.build'
    text = (S / name).read_text()
    for index in (0, 1, 2, 7):
        a, b = str(R / old['native_objects'][index]['path']), str(R / objects[index]['path'])
        assert text.count(a) == 1
        text = text.replace(a, b, 1)
    texts = {name: text}
    for name, a, b in (
        ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-inline-private-v43', 'Sirius-native-cpu-dispatch-token-v44'),
        ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8AS', 'SIR-CPU-25.2.8AT')):
        text = (S / name).read_text()
        assert text.count(a) == 1
        texts[name] = text.replace(a, b, 1)
    basis_headers={name:(S/name).read_text() for name in HEADERS}
    for name,text in basis_headers.items():
        assert data_pin(text.encode())==value['header_basis'][name]
    candidate_headers=header_patch(basis_headers,value['mesa_patch'])
    assert header_patch(candidate_headers,value['mesa_patch'],reverse=True)==basis_headers
    candidate_headers=add_layout_assertions(candidate_headers)
    assert header_patch(remove_layout_assertions(candidate_headers),value['mesa_patch'],reverse=True)==basis_headers
    texts.update(candidate_headers)
    patch, after, modified = '', dict(before), {}
    for name, text in texts.items():
        patch += ''.join(difflib.unified_diff((S / name).read_text().splitlines(True), text.splitlines(True),
            n=3, fromfile='a/' + name, tofile='b/' + name))
        after[name] = data_pin(text.encode())
        modified[name] = dict(basis=before[name], candidate=after[name])
    inputs = invariant_basis_inputs(old, value)
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
        source_inputs=inputs, new_native_compilations=0, sealed_component_native_compilations=4,
        reused_raw_objects=3, new_optional_raw_objects=0, reused_optional_raw_objects=3,
        new_component_scalar_objects=3, new_component_adapter_objects=1,
        guard_shared_byte_authority=value['guard_shared_byte_authority'],
        guard_basis_authority=old['guard_shared_byte_authority'],
        mutable_header_basis=value['header_basis'],
        admission_component_manifest=str(value['manifest'].relative_to(R)),
        admission_header_patch_identity=data_pin(value['mesa_patch'].encode()),
        observer_layout_assertion_identity=data_pin(LAYOUT_ASSERTIONS.encode()),
        header_domain_component=old['header_domain_component'],
        scope='AT replaces only baseline scalar3 and adapter; AP pure admission, oldAVX2raw3, ASraw3 and math3 remain exact. A two-header branch uses the existing positive job preflight capability for private workers that skip only repeated immutable-program comparisons and repeat all other guards. Public guarded paths, failed-preflight JIT behavior, scalar tails, FP/ISA scope,9380tuples, buffers/layouts/counts/scheduler/cache/governors remain. No lp_state_cs.c or job-layout change. Full source inverse and actual private-entry observations are required; no predecessor acceptance transfers.'))
    print(json.dumps(dict(provider_patch_files=5, new_native_compilations=0, merged_parts=7,
        reused_native_indices=[3,4,5,6], sealed_component_native_compilations=4)))
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
        scope='Reconfigure the same Mesa cache/options/toolchain for one merged seven-input object; baseline scalar3/adapter and two declared header branches changed; original pure admission/raw/math reused. No compiler/native workload or qualification claim.')
    write(P / 'actions.json', {'configure': action})
    for name in ('run_owned.py', 'owned_guard.py'):
        (P / name).write_bytes((W / name).read_bytes())
    print(json.dumps(dict(configure_inputs=len(inputs), merged_parts=7)))
else:
    raise ValueError(phase)
