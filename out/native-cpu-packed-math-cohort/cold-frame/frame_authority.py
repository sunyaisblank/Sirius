"""Pin the prepared AO frame source chain; bind AP identities and final gates."""
from pathlib import Path
import hashlib
import json

AO_PINS = {
    'frame_authority.py': '96672cde40252ad92f6871c107340281d4238c4227c5599742ebcfd19facc90a',
    'prepare.py': '95b8c9a24dcfafe92121949f7f7feae22c048f287ab4c9245685c0ac627bcadf',
    'read_result.py': '31c701f50740434be9c10724367c5cb7278046188cf79d64291b74d05d4f6b32',
    'run_owned.py': '16163059420abba7a409bebdc1b08b85c5cb85d1375af0344f6980eba62e9d5b',
    'owned_guard.py': '4d571702b673b17d5d4e95ce02ca14543bb08a47d4a14f685158985fb4d5624a',
}


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,
                    sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def ao(name):
    path = Path.cwd() / 'out/native-cpu-v4-cohort/cold-frame' / name
    data = path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == AO_PINS[name], path
    source = data.decode()
    for before, after in (
        ('native-cpu-v4-cohort', 'native-cpu-packed-math-cohort'),
        ('out/native-cpu-v4-provider/provider', 'out/native-cpu-packed-math-provider/provider'),
        ('current_AO_link_proof', 'current_AR_link_proof'),
        ('fresh AO all8/all22 controls passed', 'fresh AR all8/all22 controls passed'),
        ('AO only after all8/all22', 'AR only after all8/all22')):
        source = source.replace(before, after)
    before = ('exact AO provider with separately guarded optional v4 arithmetic and byte-preserved original v3 fallback, '
              'original AM identities before C++ calls and unchanged original general expressions/AK refinements/raw3/scalar tails/scheduler/FP scope')
    after = ('exact AR provider preserving the AP header admission extension and first7 native objects, '
             'unchanged original shader acceptance/refusal and output behavior, baseline packet selector plus same-eight-lane raw-v4 stages with private math banks and reviewed zero-residual arithmetic, '
             'separately guarded optional v4/original v3 arithmetic and original AVX2/scalar fallback/scheduler/FP scope')
    if name == 'prepare.py':
        assert source.count(before) == 1
        source = source.replace(before, after, 1)
    return source


def load(name):
    # The pinned AO authority retains and checks every original AN source pin.
    namespace = {}
    basis = Path.cwd() / 'out/native-cpu-v4-cohort/cold-frame/frame_authority.py'
    exec(compile(ao('frame_authority.py'), str(basis), 'exec'), namespace)
    source = namespace['load'](name)
    if name == 'prepare.py':
        anchor = '    inputs |= {path for path in W.iterdir() if path.is_file()}'
        assert source.count(anchor) == 1
        source = source.replace(anchor, anchor + '\n    from frame_authority import authority_files, acceptance_inputs\n'
                                '    inputs |= set(authority_files()) | acceptance_inputs()', 1)
    return source


def authority_files():
    root = Path.cwd()
    return ([root / 'out/native-cpu-v4-cohort/cold-frame' / name for name in AO_PINS] +
            [root / 'out/native-cpu-local-identities-cohort/cold-frame' / name
             for name in ('prepare.py', 'read_result.py', 'run_owned.py', 'owned_guard.py')])


def acceptance_inputs():
    # Preserve the full current sole-ABBA input and terminal-owner closure
    # in addition to the original control-7 frame prerequisites.
    root = Path.cwd()
    benefit = root / 'out/native-cpu-packed-math-cohort/benefit'
    actions = benefit / 'runtime/actions.json'
    cohort = json.loads(actions.read_text())
    inputs = {root / name for name in cohort['b2']['frozen_inputs']}
    inputs.add(actions)
    for label in ('a1', 'b1', 'b2', 'a2'):
        for name in ('owner.json', 'gtest.xml'):
            inputs.add(benefit / 'runtime' / label / name)
    return inputs


def check_acceptance():
    """Do not prepare or start a frame before the current AR final admission."""
    root = Path.cwd()
    current = root / 'out/native-cpu-packed-math-cohort'
    paths = dict(controls=current / 'controls/readers/controls-1-8-acceptance.json',
                 benefit=current / 'benefit/readers/benefit.json',
                 decimal=current / 'benefit/readers/decimal-agreement.json',
                 c8=current / 'control-8-acceptance.json')
    rows = {name: json.loads(path.read_text()) for name, path in paths.items()}
    controls, benefit, decimal, c8 = (rows[name] for name in ('controls', 'benefit', 'decimal', 'c8'))
    assert controls['all8_original_controls_pass'] and controls['all22_original_benefit_checks_pass']
    assert benefit['audit_passed'] and benefit['all22_numeric_pass'] and benefit['checks_passed'] == 22
    assert decimal['agrees_on_all22'] and decimal['all22_pass']
    assert c8['passed'] and c8['all32_readbacks_match_preserved_originals']
    assert c8['optional_math']['all_selected_optional_math_observed']
    for key, name in (('control8_acceptance', 'c8'), ('benefit', 'benefit'), ('decimal', 'decimal')):
        assert controls[key] == pin(paths[name]), key
    for label in ('a1', 'b1', 'b2', 'a2'):
        owner = json.loads((current / 'benefit/runtime' / label / 'owner.json').read_text())
        assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
