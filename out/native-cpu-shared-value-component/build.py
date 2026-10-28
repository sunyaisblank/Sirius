"""Compile the admitted adapter once and bind the shared core object once."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-shared-value-component'
V = R / 'out/native-cpu-strict-ir-linking'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,
                    sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


phase, = sys.argv[1:]
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    derivation = json.loads((W / 'derivation.json').read_text())
    assert pin(W / 'rows_adapter.cpp') == {k: derivation['derived'][k] for k in ('bytes', 'sha256')}
    previous = json.loads((R / 'out/native-cpu-dispatch-token-provider/merge-inputs.json').read_text())
    reused = previous['parts'][1:]
    assert len(reused) == 6
    parts = []
    for stage in ('dp', 'transport', 'endpoint', 'math'):
        case = 'shared-endpoint-cleanup-native' if stage == 'endpoint' else 'shared-' + stage + '-native'
        owner = json.loads((V / case / 'owner.json').read_text())
        assert owner['passed'] and owner['status'] == 'terminal' and owner['whole_inputs_unchanged']
        product = json.loads((V / (case + '-result.json')).read_text())['output']
        parts.append(product)
    for row in reused + parts:
        assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}
    compile_argv = list(derivation['compile_argv'])
    compile_argv[-1] = str(B / 'adapter-only.o')
    merge_argv = ['/usr/bin/ld', '-r', '-o', str(B / 'rows_adapter.o'), str(B / 'adapter-only.o')]
    merge_argv += [str(R / row['path']) for row in reused + parts]
    sources = [*W.glob('*.py'), W / 'rows_adapter.cpp', W / 'derivation.json',
        R / 'out/native-cpu-dispatch-token-provider/merge-inputs.json',
        R / 'out/native-cpu-dispatch-admission-token/admitted_bridge.h',
        R / 'out/native-cpu-dispatch-admission-token/native_packet.hpp',
        R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
        R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py']
    sources += [R / row['path'] for row in reused + parts]
    sources += [Path(compile_argv[0]).resolve(), Path('/usr/bin/ld').resolve()]
    inputs = {str(p.relative_to(R)) if p.is_relative_to(R) else str(p): pin(p) for p in sources}
    write(W / 'build-plan.json', dict(inputs=inputs, compile_argv=compile_argv, merge_argv=merge_argv,
        reused=reused, new_value_objects=parts,
        scope='One new baseline adapter compile plus11-part relocation link. Public wrappers and6 old AS/math/dispatcher objects retained, only private admitted v4 branches invoke new scoped raw3/shared9. No numerical/runtime claim.'))
    assert not B.exists()
    write(W / 'actions.json', {'build': dict(
        revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=['/usr/bin/python3', '-B', str(W / 'build.py'), 'build'],
        frozen_inputs=sorted(inputs) + ['out/native-cpu-shared-value-component/build-plan.json'],
        scope='One baseline adapter compile and exact11-object merged native link; original public branches, scalar/AVX2/oldv4 math retained. Static relocation/ABI evidence only; no numerical, provider, frame or release acceptance.')})
elif phase == 'build':
    plan = json.loads((W / 'build-plan.json').read_text())
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    B.mkdir()
    subprocess.run(plan['compile_argv'], check=True)
    subprocess.run(plan['merge_argv'], check=True)
    target = B / 'rows_adapter.o'
    undefined = subprocess.check_output(['nm', '-u', str(target)], text=True)
    assert 'sirius_shared_value_' not in undefined and 'sirius_scope_' not in undefined, undefined
    listing = subprocess.check_output(['readelf', '-Ws', str(target)], text=True)
    cores = []
    for line in listing.splitlines():
        fields = line.split()
        if len(fields) == 8 and fields[-1].startswith('sirius_shared_value_'):
            assert fields[3:6] == ['FUNC', 'GLOBAL', 'HIDDEN'] and fields[6] != 'UND'
            assert int(fields[2]) > 0
            cores.append(fields[-1])
    expected = json.loads((V / 'shared-plan.json').read_text())['core_names']
    assert sorted(cores) == sorted(expected)
    sources = []
    for row in plan['new_value_objects']:
        refs = subprocess.check_output(['nm', '-u', str(R / row['path'])], text=True)
        sources.append(dict(object=row, shared_core_references=[line.split()[-1] for line in refs.splitlines() if 'sirius_shared_value_' in line]))
    write(W / 'build-result.json', dict(passed=True,
        output=dict(path=str(target.relative_to(R)), **pin(target)), adapter=pin(B / 'adapter-only.o'),
        shared_core_definitions=cores, all_raw_scope_and_core_references_bound=True,
        raw_object_references=sources, remaining_external_symbols=undefined,
        reused=plan['reused'], scope='Exact merged object/static symbol joins only. Full provider PC/ISA/FP-scope, whole numerical, C8/all32/ABBA/all8 and original frame qualification remain required.'))
    print(json.dumps(dict(passed=True, shared_core_definitions=len(cores), merged_bytes=target.stat().st_size)))
else:
    raise ValueError(phase)
