"""Emit the two remaining original raw modules once, using the same compiler."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-strict-ir-linking'


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
    proof = json.loads((R / 'out/native-cpu-inline-private-bank/derivation.json').read_text())
    source_inputs = dict(proof['source_inputs']) | dict(proof['derived'])
    for name, identity in source_inputs.items():
        assert pin(R / name) == identity, name
    for path in (*W.glob('*.py'), R / 'out/native-cpu-inline-private-bank/derivation.json',
                 R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
                 R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py'):
        source_inputs[str(path.relative_to(R))] = pin(path)
    actions = json.loads((W / 'actions.json').read_text())
    for stage in ('transport', 'endpoint'):
        name = 'emit-' + stage
        assert name not in actions and not (B / (stage + '.ll')).exists()
        original, = [row for row in proof['raw_compilations'] if row['stage'] == stage]
        argv = list(original['argv'])
        assert argv[-3:-1] == ['--emit-obj', '-o']
        argv[-3], argv[-1] = '--emit-llvm-text', str(B / (stage + '.ll'))
        plan_path = W / (stage + '-plan.json')
        write(plan_path, dict(argv=argv, source_inputs=source_inputs,
                              only_output_kind_and_path_changed=True))
        actions[name] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            kind='diagnostic', expected_exit=0, software_provider=False,
            limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
            argv=['/usr/bin/python3', '-B', str(W / 'emit_remaining.py'), stage],
            frozen_inputs=sorted(source_inputs) + [str(plan_path.relative_to(R))],
            scope='First LLVM emission of this unchanged original AS stage with identical ISPC/O2/SKX/PIC/disable-fma/system-math flags. Required vocabulary input for coupled3-stage implementation; no conversion, runtime, numerical, speed or adoption claim.')
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=['emit-transport', 'emit-endpoint'], compiler_run=False)))
else:
    assert phase in ('transport', 'endpoint')
    # Reuse the exact already-executed inventory implementation, not a new
    # interpretation. Adapt only its plan/output/receipt paths and existing-B
    # assertion. This reuses code, not DP evidence for a different module.
    basis = W / 'emit_ir.py'
    source = basis.read_text()
    assert hashlib.sha256(source.encode()).hexdigest() == '0f09f0f925db44daaec7dd4dc8c151a4d816627aff13fa9e5e66b3b6b1c4172b'
    prefix = source[:source.index('phase, = sys.argv[1:]')]
    body = source.split("elif phase == 'emit':\n", 1)[1].rsplit('\nelse:\n', 1)[0]
    import textwrap
    body = textwrap.dedent(body)
    replacements = {"'plan.json'": repr(phase + '-plan.json'),
                    "B.mkdir(exist_ok=False)": 'assert B.is_dir()',
                    "'dopri_phase.ll'": repr(phase + '.ll'),
                    "'inventory.json'": repr(phase + '-inventory.json')}
    for old, new in replacements.items():
        assert body.count(old) == 1, old
        body = body.replace(old, new, 1)
    scope = {'__file__': str(basis), '__name__': 'remaining_raw_ir_inventory'}
    exec(compile(prefix + body, str(basis), 'exec'), scope)
