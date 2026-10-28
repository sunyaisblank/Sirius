"""Prepare/run one12-input relocatable merge; compile no adapter or raw source."""
from pathlib import Path
import json
import subprocess
import sys
from component import repaired, doc, pin, checked, CORES, LOCAL_CORES

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-closed-dp-unwind-component'


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


phase, = sys.argv[1:]
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    value = repaired()
    assert not B.exists() and not (W / 'actions.json').exists()
    parts = value['parts']
    argv = ['/usr/bin/ld', '-r', '-o', str(B / 'rows_adapter.o')] + [str(R / row['path']) for row in parts]
    inputs = dict(value['inputs'])
    for path in [*W.glob('*.py'), Path('/usr/bin/ld').resolve(), Path('/usr/bin/nm').resolve(), Path('/usr/bin/readelf').resolve()]:
        inputs[str(path.relative_to(R)) if path.is_relative_to(R) else str(path)] = pin(path)
    write(W / 'build-plan.json', dict(inputs=inputs, parts=parts, merge_argv=argv, new_native_compilations=0,
        adapter_recompiled=False, dp_unwind_repair=value['unwind_repair'],
        scope='One ld-r only; reuse eleven AV constituents including the adapter and replace only DP with its code-identical unwind-table repair. Static interface only, no numerical transfer.'))
    write(W / 'actions.json', {'build': dict(
        revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=['/usr/bin/python3', '-B', str(W / 'build.py'), 'build'],
        frozen_inputs=sorted(inputs) + ['out/native-cpu-closed-dp-unwind-component/build-plan.json'],
        scope='One12-input relocatable link, no GCC/ISPC/LLVM compilation, provider loading or workload. New FDE metadata only; current C8 acceptance is not inherited.')})
elif phase == 'build':
    value = repaired()
    plan = doc(W / 'build-plan.json')
    assert plan['parts'] == value['parts'] and plan['new_native_compilations'] == 0
    for name, expected in plan['inputs'].items():
        checked(R / name, expected)
    assert not B.exists()
    B.mkdir()
    subprocess.run(plan['merge_argv'], check=True)
    target = B / 'rows_adapter.o'
    undefined = subprocess.check_output(['/usr/bin/nm', '-u', str(target)], text=True)
    assert not any(prefix in undefined for prefix in ('sirius_shared_value_', 'sirius_local_value_', 'sirius_scope_', 'sirius_closed_'))
    listing = subprocess.check_output(['/usr/bin/readelf', '-Ws', str(target)], text=True)
    wanted = CORES | LOCAL_CORES
    cores, closed = [], []
    for line in listing.splitlines():
        fields = line.split()
        if len(fields) != 8 or fields[-1] not in wanted | {'sirius_closed_dp_packet'}:
            continue
        size = int(fields[2], 16 if fields[2].startswith('0x') else 10)
        assert fields[3:6] == ['FUNC', 'GLOBAL', 'HIDDEN'] and fields[6] != 'UND' and size > 0
        (cores if fields[-1] in wanted else closed).append(fields[-1])
    assert set(cores) == wanted and len(cores) == 18 and closed == ['sirius_closed_dp_packet']
    write(W / 'build-result.json', dict(passed=True, output=dict(path=str(target.relative_to(R)), **pin(target)),
        parts=value['parts'], adapter=value['parts'][0], core_definitions=sorted(cores), closed_raw_exports=closed,
        all_raw_scope_and_core_references_bound=True, compiled_interface_passed=True, new_native_compilations=0,
        numerical_qualification_transferred=False, remaining_external_symbols=undefined,
        scope='One ld-r/static symbol gate only; no adapter/source/math rebuild and no inherited numerical or C8 verdict. Final DSO code/FDE and genuine-job original C8 remain required.'))
else:
    raise ValueError(phase)
