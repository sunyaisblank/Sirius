"""One fixed dyadic/boundary matrix against the sealed original arithmetic."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-exact-product-classes'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def write(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')


def terminal(name):
    value = doc(W / name / 'owner.json')
    assert value['passed'] and value['status'] == 'terminal' and value['child_exit'] == 0
    assert value['source_unchanged'] and value['whole_inputs_unchanged']
    assert not value['cleanup_errors'] and not value['remaining_owned_processes']


phase, = sys.argv[1:]
proof = doc(W / 'derivation.json')
assert doc(W / 'compatibility.json')['passed']
assert pin(B / 'optional.o') == doc(W / 'compile-result.json')['optional_object']
for row in proof['reused_objects']:
    assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}
if phase in ('prepare-build', 'prepare-run'):
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    actions = doc(W / 'actions.json')
    inputs = set(actions['run']['frozen_inputs'])
    inputs |= {str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
    if phase == 'prepare-build':
        args = list(proof['compile_argv'])
        args = args[:args.index('-c')]
        args[args.index('-march=x86-64-v4')] = '-march=x86-64'
        args.remove('-mprefer-vector-width=256')
        args += ['-c', str(W / 'driver.cpp'), '-o', str(B / 'matrix.o')]
        link = [args[0], str(B / 'matrix.o'),
                str(R / 'bin/linux-gcc/native-cpu-v4-arithmetic/dispatch.o'),
                str(R / 'bin/linux-gcc/native-cpu-v4-arithmetic/v3.o'),
                str(B / 'optional.o'), '-Wl,--no-relax', '-o', str(B / 'matrix-check')]
        write(W / 'matrix-plan.json', dict(compile=args, link=link))
        inputs.add(str((W / 'matrix-plan.json').relative_to(R)))
        name, argv = 'matrix-build', ['/usr/bin/python3', '-B', str(W / 'matrix.py'), 'build']
        scope = 'One baseline fixed-case driver compile/link, sealed old v3/dispatcher and new optional object. Integer-only constexpr center identities; no scalar reference or raw/provider recompilation.'
    else:
        terminal('matrix-build')
        inputs |= {str((B / n).relative_to(R)) for n in ('matrix.o', 'matrix-check')}
        name, argv = 'matrix-run', [str(B / 'matrix-check')]
        scope = 'One fixed17-packet/136-lane matrix,17actual dispatcher plus17sealed-v3 calls.2040bank/272canary words, immutable inputs, nine independent dyadic center cases, strict MXCSR flags/state and actual denormal reference witness. No random cases, benchmark, branch-PC, hit-rate, provider/frame/scientific/release admission.'
    assert name not in actions and not (W / name).exists()
    actions[name] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=argv, frozen_inputs=sorted(inputs), scope=scope)
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=name, executed=False)))
elif phase == 'build':
    for key in ('compile', 'link'):
        subprocess.run(doc(W / 'matrix-plan.json')[key], check=True)
elif phase == 'read':
    terminal('matrix-run')
    result = json.loads((W / 'matrix-run/stdout.log').read_text())
    assert result['passed'] and all(row['passed'] for row in result['cases'])
    for key, expected in dict(case_count=17, candidate_calls=17, original_v3_reference_calls=17,
            active_lane_results=136, compared_bank_words=2040, compared_canary_words=272,
            independent_center_cases=9, independent_center_lane_results=72,
            reference_only_cases=8, immutable_input_words_per_case=80).items():
        assert result[key] == expected, key
    assert sum(row['independent_center_checked'] for row in result['cases']) == 9
    assert result['subnormal_power_reference_denormal_operand_observed']
    write(W / 'matrix-compatibility.json', dict(passed=True, comparison=result,
        optional_object=pin(B / 'optional.o'), executable=pin(B / 'matrix-check'),
        scope='Finite fixed-case exact reference/center/FP-state compatibility only; no measured performance, classifier PC coverage, provider/frame or scientific/release acceptance.'))
    print(json.dumps(dict(passed=True, cases=17, bank_words=2040,
                          independent_centers=72, exact_MXCSR_flags_and_state=True)))
else:
    raise ValueError(phase)
