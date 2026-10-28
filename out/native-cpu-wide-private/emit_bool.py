"""Correct ISPC bool-select syntax, retaining the first failed emit and source owner."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
RUN = W / 'raw-bool'
BASE = W / 'emit_raw.py'
BASE_SHA = 'b368905160fcceb16ea6e803507fa2c993f2c491acd675af54265717e7a9854e'
OLD_REV = '258f1dda0849aed95fbcbaa6d15c74bca693cea3'
HELPER = W / 'subgroup8.isph'
ORDER = ('emit-transport', 'emit-endpoint', 'emit-dopri_phase')


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


assert pin(BASE)['sha256'] == BASE_SHA
phase, = sys.argv[1:]
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    assert not (RUN / 'actions.json').exists()
    old = subprocess.check_output(['git', 'show', OLD_REV + ':out/native-cpu-wide-private/subgroup8.isph'])
    before = doc(W / 'raw/emit-transport/inputs-before.json')
    key = str(HELPER.relative_to(R))
    assert hashlib.sha256(old).hexdigest() == before[key]['sha256']
    expected = old.decode()
    edits = [('select(active, value, true)', 'active ? value : true', 1),
             ('select(active, value, false)', 'active ? value : false', 1),
             ('select(programIndex < 8, low, high)', '(programIndex < 8) ? low : high', 2)]
    for original, corrected, count in edits:
        assert expected.count(original) == count and corrected not in expected
        expected = expected.replace(original, corrected)
    assert HELPER.read_text() == expected
    failed = doc(W / 'raw/emit-transport/owner.json')
    assert failed['status'] == 'terminal' and not failed['passed'] and failed['child_exit'] == 1
    assert failed['source_unchanged'] and failed['whole_inputs_unchanged']
    assert not failed['remaining_owned_processes'] and not failed['cleanup_errors']
    assert failed['observed_births'] and all(row['absent'] for row in failed['observed_births'])
    assert (W / 'raw/emit-transport/inputs-before.json').read_bytes() == (W / 'raw/emit-transport/inputs-after.json').read_bytes()
    diagnostic = (W / 'raw/emit-transport/stderr.log').read_text()
    assert 'Multiple overloaded functions matched call to function' in diagnostic
    assert not any((W / 'raw' / name).exists() for name in ORDER[1:])
    plan = doc(W / 'raw-plan.json')
    inputs = dict(plan['inputs'])
    assert inputs[key] == {k:before[key][k] for k in ('bytes', 'sha256')}
    inputs[key] = pin(HELPER)
    source = doc(W / 'raw-source-contract.json')
    assert source['outputs'] == plan['expected_source_outputs']
    for name, expected in source['outputs'].items():
        assert pin(R / name) == expected
        inputs[name] = expected
    extras = [BASE, Path(__file__).resolve(), RUN / 'run_owned.py', RUN / 'owned_guard.py',
              W / 'raw-source-contract.json', W / 'derive-result.json', W / 'raw/derive/owner.json',
              W / 'raw-plan.json', W / 'raw/actions.json']
    extras += [W / 'raw/emit-transport' / name for name in ('owner.json', 'inputs-before.json', 'inputs-after.json', 'stdout.log', 'stderr.log')]
    for path in extras:
        inputs[str(path.relative_to(R))] = pin(path)
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    commands = {name:plan['commands'][name] for name in ORDER}
    assert all(not (R / argv[-1]).exists() for argv in commands.values())
    correction = dict(original_revision=OLD_REV, original_source=old.decode(),
        original_pin={k:before[key][k] for k in ('bytes', 'sha256')}, corrected_pin=pin(HELPER),
        exact_replacements=[dict(before=a, after=b, count=c) for a,b,c in edits],
        failed_owner=pin(W / 'raw/emit-transport/owner.json'),
        source_derivation_reused=True, numerical_or_mask_semantics_changed=False)
    write(W / 'raw-bool-plan.json', dict(inputs=inputs, commands=commands,
        expected_source_outputs=source['outputs'], math_interface=plan['math_interface'], correction=correction))
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    actions = {}
    for i, name in enumerate(ORDER):
        previous = [] if not i else [str((RUN / ORDER[i-1] / 'owner.json').relative_to(R)),
            str((W / (ORDER[i-1] + '-result.json')).relative_to(R)),
            str((Path(commands[ORDER[i-1]][-1])).relative_to(R))]
        actions[name] = dict(revision=revision, kind='diagnostic', expected_exit=0, software_provider=False,
            limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
            argv=['/usr/bin/python3', '-B', str(Path(__file__).resolve()), name],
            frozen_inputs=sorted(inputs) + [str((W / 'raw-bool-plan.json').relative_to(R))] + previous,
            scope='First corrected raw16 emit after one preserved bool-select API failure. Reuse the exact '
                  'successful source derivation and strict8 adapter objects; four bool conditional expressions '
                  'replace unsupported overloaded select calls. No changed arithmetic, masks, scientific work '
                  'or limits, no repeated program study, and no strict/native/runtime acceptance.')
    write(RUN / 'actions.json', actions)
    print(json.dumps(dict(prepared=list(ORDER), reused_source_derivation=True, preserved_failed_emit=True)))
else:
    assert phase in ORDER
    source = BASE.read_text()
    changes = [("RUN = W / 'raw'", "RUN = W / 'raw-bool'"),
               ("W / 'raw-plan.json'", "W / 'raw-bool-plan.json'"),
               ("terminal(RUN / 'derive')", "terminal(W / 'raw/derive')")]
    for original, corrected in changes:
        assert source.count(original) >= 1
        source = source.replace(original, corrected)
    anchor = "    write(W / (phase + '-result.json'), result)"
    assert source.count(anchor) == 1
    source = source.replace(anchor,
        "    result['accepted_artifact_owner'] = str((RUN / phase / 'owner.json').relative_to(R))\n"
        "    result['source_correction'] = pin(W / 'raw-bool-plan.json')\n" + anchor, 1)
    exec(compile(source, str(BASE), 'exec'))
