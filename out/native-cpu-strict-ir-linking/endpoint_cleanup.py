"""Apply only the required interface cleanup to already-O2 Endpoint IR."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import textwrap

basis = Path.cwd() / 'out/native-cpu-strict-ir-linking/shared.py'
source = basis.read_text()
assert hashlib.sha256(source.encode()).hexdigest() == '09aa472bb8916df5bf6ab4a469f148d61f6fb375031e7c38bb383358763a1e9b'
scope = {'__file__': str(basis), '__name__': 'focused_endpoint_cleanup'}
exec(compile(source[:source.index('phase, = sys.argv[1:]')], str(basis), 'exec'), scope)
R, W, B, D, pin, write = (scope[k] for k in ('R', 'W', 'B', 'D', 'pin', 'write'))
phase, = sys.argv[1:]
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    failed = json.loads((W / 'shared-endpoint-opt/owner.json').read_text())
    assert failed['status'] == 'terminal' and failed['stop_reason'] == 'diagnostic_timeout'
    assert not failed['remaining_owned_processes'] and not failed['cleanup_errors']
    lower = json.loads((W / 'shared-endpoint-lower/owner.json').read_text())
    assert lower['passed'] and lower['status'] == 'terminal'
    product = json.loads((W / 'shared-endpoint-lower-result.json').read_text())['output']
    assert pin(R / product['path']) == {k: product[k] for k in ('bytes', 'sha256')}
    outputs = {'shared-endpoint-cleanup-opt': 'shared-endpoint.cleanup.O2.bc',
               'shared-endpoint-cleanup-native': 'shared-endpoint.cleanup.o'}
    pipeline = 'always-inline,function(sroa,instcombine,simplifycfg,early-cse<memssa>),globaldce,verify'
    commands = {
        'shared-endpoint-cleanup-opt': [str(D / 'opt'), '-passes=' + pipeline, '-verify-each',
            '-mcpu=skylake-avx512', str(R / product['path']), '-o', str(B / outputs['shared-endpoint-cleanup-opt'])],
        'shared-endpoint-cleanup-native': [str(D / 'llc'), '-O2', '-filetype=obj', '-relocation-model=pic',
            '-mcpu=skylake-avx512', str(B / outputs['shared-endpoint-cleanup-opt']),
            '-o', str(B / outputs['shared-endpoint-cleanup-native'])]}
    inputs = dict(json.loads((W / 'shared-plan.json').read_text())['inputs'])
    for p in (Path(__file__).resolve(), W / 'shared-endpoint-opt/owner.json',
              W / 'shared-endpoint-opt/stderr.log', W / 'shared-endpoint-lower-result.json',
              R / product['path']):
        inputs[str(p.relative_to(R))] = pin(p)
    write(W / 'endpoint-cleanup-plan.json', dict(commands=commands, outputs=outputs, inputs=inputs,
        reason='Original AS input is already ISPC-O2. The new change needs always-inline adapters and SROA of initialized private banks, followed by local cleanup; avoid a second full whole-module O2 pipeline.',
        prior_timeout_preserved=True, work_limits_unchanged=True))
    actions = json.loads((W / 'actions.json').read_text())
    for name in commands:
        assert name not in actions and not (B / outputs[name]).exists()
        actions[name] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            kind='diagnostic', expected_exit=0, software_provider=False,
            limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
            argv=['/usr/bin/python3', '-B', str(Path(__file__).resolve()), name],
            frozen_inputs=sorted(inputs) + ['out/native-cpu-strict-ir-linking/endpoint-cleanup-plan.json'] +
                ([str((B / outputs['shared-endpoint-cleanup-opt']).relative_to(R))] if name.endswith('-native') else []),
            scope='Focused cleanup of the changed interface in already-O2 Endpoint IR; same structural/signature/strict-FP/no-FMA gates and unchanged work limits. No numerical or performance/provider/frame acceptance.')
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=list(commands), pipeline=pipeline, compiler_run=False)))
else:
    assert phase in ('shared-endpoint-cleanup-opt', 'shared-endpoint-cleanup-native')
    anchor = "else:\n    plan = json.loads((W / 'shared-plan.json').read_text())"
    assert source.count(anchor) == 1
    body = textwrap.dedent(source.split(anchor, 1)[1])
    body = "plan = json.loads((W / 'endpoint-cleanup-plan.json').read_text())\n" + body
    old = "dependency = 'shared-plugin' if phase.endswith('-lower') else phase.rsplit('-', 1)[0] + ('-lower' if phase.endswith('-opt') else '-opt')"
    assert body.count(old) == 1
    body = body.replace(old, "dependency = 'shared-endpoint-lower' if phase.endswith('-opt') else 'shared-endpoint-cleanup-opt'", 1)
    scope['phase'] = phase
    exec(compile(body, str(basis), 'exec'), scope)
