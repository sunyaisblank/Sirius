"""Own one positive fixture build and one passive fixed-case observation."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-exact-product-residual'


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
compiled = doc(W / 'compile-result.json')
assert pin(B / 'optional.o') == compiled['optional_object']
for row in proof['reused_objects']:
    assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}
if phase in ('prepare-build', 'prepare-run'):
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    actions = doc(W / 'actions.json')
    inputs = set(actions['run']['frozen_inputs'])
    inputs |= {str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
    inputs |= {'/usr/bin/gdb', 'out/native-cpu-accumulate-compatibility/elf_proof.py'}
    if phase == 'prepare-build':
        args = list(proof['compile_argv'])
        args = args[:args.index('-c')]
        args[args.index('-march=x86-64-v4')] = '-march=x86-64'
        args.remove('-mprefer-vector-width=256')
        args += ['-c', str(W / 'positive.cpp'), '-o', str(B / 'positive.o')]
        link = [args[0], str(B / 'positive.o'),
                str(R / 'bin/linux-gcc/native-cpu-v4-arithmetic/dispatch.o'),
                str(R / 'bin/linux-gcc/native-cpu-v4-arithmetic/v3.o'),
                str(B / 'optional.o'), '-Wl,--no-relax', '-o', str(B / 'positive-check')]
        write(W / 'positive-plan.json', dict(compile=args, link=link))
        inputs.add(str((W / 'positive-plan.json').relative_to(R)))
        name = 'positive-build'
        argv = ['/usr/bin/python3', '-B', str(W / 'positive_owned.py'), 'build']
        scope = 'One baseline positive-fixture compile and exact old-v3/dispatch/new-optional link. One fixed eight-lane dyadic product, not a benchmark or shader/provider trial.'
    else:
        terminal('positive-build')
        inputs |= {str((B / n).relative_to(R)) for n in ('positive.o', 'positive-check')}
        name = 'positive-run'
        argv = ['/usr/bin/gdb', '-q', '-nx', '-batch', '-x', str(W / 'positive.gdb'), str(B / 'positive-check')]
        scope = 'One fixed eight-lane positive product via real dispatcher versus sealed v3, complete bank/canaries/immutable inputs, independent dyadic center and exact MXCSR flags/state. Passive first-residual PC/register read only; no inferior calls, data/register writes, retries, sampling for success or timing/adoption claim.'
    assert name not in actions and not (W / name).exists()
    actions[name] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=argv, frozen_inputs=sorted(inputs), scope=scope)
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=name, executed=False)))
elif phase == 'build':
    for key in ('compile', 'link'):
        subprocess.run(doc(W / 'positive-plan.json')[key], check=True)
elif phase == 'read':
    terminal('positive-run')
    observed = doc(W / 'positive-run/targets.json')
    assert observed['coverage_passed'] and observed['exit_code'] == 0
    result, = [json.loads(line) for line in (W / 'positive-run/stdout.log').read_text().splitlines()
               if line.startswith('{"passed"')]
    assert result['passed'] and result['candidate_calls'] == result['original_v3_reference_calls'] == 1
    assert result['active_lane_results'] == 8 and result['compared_bank_words'] == 120
    assert result['MXCSR_sticky_equal'] and result['caller_MXCSR_state_restored']
    write(W / 'positive-compatibility.json', dict(passed=True, comparison=result,
        actual_residual_observation=observed,
        scope='One actual fast-branch positive case plus original finite corpus only. Hit rate, throughput, older-hardware execution, provider/frame/scientific/release acceptance remain unqualified.'))
    print(json.dumps(dict(passed=True, actual_fast_branch_observed=True, compared_bank_words=120,
                          private_MXCSR_sticky_equal=True)))
else:
    raise ValueError(phase)
