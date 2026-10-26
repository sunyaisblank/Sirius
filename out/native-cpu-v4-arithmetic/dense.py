"""Reuse the unchanged five-op dense fixture; passively observe first targets."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-v4-arithmetic'
OLD = R / 'out/native-cpu-direct-compound'


def pin(p):
    with p.open('rb') as f:
        return dict(bytes=p.stat().st_size, sha256=hashlib.file_digest(f, 'sha256').hexdigest())


phase, = sys.argv[1:]
prepared = json.loads((OLD / 'prepared-cpp.json').read_text())
built = json.loads((OLD / 'cpp-build-result.json').read_text())
assert built['passed']
for n, expected in prepared['inputs'].items():
    assert pin(R / n) == expected, n
driver = next(row['object'] for row in built['compilations'] if row['name'] == 'driver')
assert pin(R / driver['path']) == {k: driver[k] for k in ('bytes', 'sha256')}
if phase == 'read':
    owner = json.loads((W / 'dense-run/owner.json').read_text())
    assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    observation = json.loads((W / 'dense-run/targets.json').read_text())
    assert observation['exit_code'] == 0 and not observation['v3_entries']
    assert len(observation['v4_entries']) == 5
    result, = [json.loads(line) for line in (W / 'dense-run/stdout.log').read_text().splitlines() if line.startswith('{"passed"')]
    assert result == dict(passed=True, calls=185220, compared_words=22226400,
                         active_lane_results=555660, case_tuples=49,
                         input_alias_and_inactive_words_checked=True, FP_controls_unchanged=True)
    out = W / 'dense-compatibility.json'
    with out.open('x') as f:
        json.dump(dict(result=result, observation=observation,
             scope='Original dense scalar compatibility and five actual optional entry PCs only; '
                   'no non-v4-host execution, raw stage, benefit, frame or release qualification.'), f, indent=2)
        f.write('\n')
    print(json.dumps(result))
    raise SystemExit(0)

assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
inputs = set(prepared['inputs']) | {str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
inputs |= {driver['path']} | {str((B / n).relative_to(R)) for n in ('optional.o', 'dispatch.o', 'v3.o')}
if phase == 'link':
    argv = [prepared['link'][0], str(R / driver['path']), str(B / 'dispatch.o'),
            str(B / 'v3.o'), str(B / 'optional.o'), '-o', str(B / 'dense-check')]
    scope = 'One link of unchanged original five-op dense scalar driver; no driver recompilation.'
elif phase == 'run':
    owner = json.loads((W / 'dense-link/owner.json').read_text())
    assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    inputs |= {str((B / 'dense-check').relative_to(R)), '/usr/bin/gdb'}
    argv = ['/usr/bin/gdb', '-q', '-nx', '-batch', '-x', str(W / 'targets.gdb'), str(B / 'dense-check')]
    scope = 'One unchanged185220-call five-op dense fixture with first-entry-only passive breakpoints. '
    scope += 'Original bank/mask/alias/refusal/FP-control assertions; no inferior function calls, memory edits or timing claim.'
else:
    raise ValueError(phase)
actions = json.loads((W / 'actions.json').read_text())
name = 'dense-' + phase
assert name not in actions
actions[name] = dict(revision=revision, kind='diagnostic', expected_exit=0,
    limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152), argv=argv,
    frozen_inputs=sorted(inputs), software_provider=False, scope=scope)
(W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps(dict(prepared=name, driver_recompiled=False)))
