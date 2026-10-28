"""Pin the prepared AO frame source chain; bind AP identities and final gates."""
from pathlib import Path
import hashlib

_AS_SOURCE = Path.cwd() / 'out/native-cpu-inline-private-cohort/cold-frame'
_AS_PINS = {'frame_authority.py': 'd77bb7feea8eccdd1728e6ebbc2ade2946d8ccc6242c8a81594b4a27dce026a5', 'owned_guard.py': '5c92a9c3d087e16bed596d74e240d557994e87d9bd1351bc4206ecc0c757214c', 'prepare.py': 'ef8eae3716c55c0f9b29494cddcda59c38f08b30751be7500587d39f6d89643e', 'read_result.py': 'b71c168987ab5d34bf795c8f1447675ec8d359e7180573d6f83499440c804a86', 'run_owned.py': 'f8456383160f96e83cb5d56900b74822c09a154f36893aa28865001a0db1bbfd'}
for _name, _digest in _AS_PINS.items():
    assert hashlib.sha256((_AS_SOURCE / _name).read_bytes()).hexdigest() == _digest, _name

_AR_SOURCE = Path.cwd() / 'out/native-cpu-packed-math-cohort/cold-frame'
_AR_PINS = {'frame_authority.py': 'f8b25c4a4399793c471525c990a674540348b8833ce9eb92e73c132a76cf474b', 'owned_guard.py': '5c92a9c3d087e16bed596d74e240d557994e87d9bd1351bc4206ecc0c757214c', 'prepare.py': '1cd5106ffd2779142ad1849438bbeceb43518cf8f57c85974b60da8f75923980', 'read_result.py': 'b71c168987ab5d34bf795c8f1447675ec8d359e7180573d6f83499440c804a86', 'run_owned.py': 'cf303693aeef445e8d5a3b1c671d1473795b612fe9a97d3522030e6d0ecd2354'}
for _name, _digest in _AR_PINS.items():
    assert hashlib.sha256((_AR_SOURCE / _name).read_bytes()).hexdigest() == _digest, _name
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
        ('native-cpu-v4-cohort', 'native-cpu-dispatch-token-cohort'),
        ('out/native-cpu-v4-provider/provider', 'out/native-cpu-shared-value-provider/provider'),
        ('current_AO_link_proof', 'current_AU_link_proof'),
        ('fresh AO all8/all22 controls passed', 'fresh AU all8/all22 controls passed'),
        ('AO only after all8/all22', 'AT only after all8/all22')):
        source = source.replace(before, after)
    before = ('exact AO provider with separately guarded optional v4 arithmetic and byte-preserved original v3 fallback, '
              'original AM identities before C++ calls and unchanged original general expressions/AK refinements/raw3/scalar tails/scheduler/FP scope')
    after = ('exact AU provider retaining AT public/scalar/pure admission and all AS public raw/GCC math/AVX2/v3 fallbacks; '
             'only genuine positive-job v4 private callbacks enter new strict raw3 and shared-nine-core value math inside the unchanged full saved/masked/restored MXCSR scope, '
             'same original actual-Vulkan host/scene/counts/layouts/scheduler/resources and all original numerical, control and work limits')
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


_previous_authority_files = authority_files
def authority_files():
    return _previous_authority_files() + [_AR_SOURCE / name for name in _AR_PINS]


def acceptance_inputs():
    # Preserve the full current sole-ABBA input and terminal-owner closure
    # in addition to the original control-7 frame prerequisites.
    root = Path.cwd()
    benefit = root / 'out/native-cpu-shared-value-cohort/benefit'
    actions = benefit / 'runtime/actions.json'
    cohort = json.loads(actions.read_text())
    inputs = {root / name for name in cohort['b2']['frozen_inputs']}
    inputs.add(actions)
    for label in ('a1', 'b1', 'b2', 'a2'):
        for name in ('owner.json', 'gtest.xml'):
            inputs.add(benefit / 'runtime' / label / name)
    return inputs


def check_acceptance():
    """Do not prepare or start a frame before the current AU final admission."""
    root = Path.cwd()
    current = root / 'out/native-cpu-shared-value-cohort'
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
    assert c8['optional_math']['genuine_private_FP_scope_observed']
    assert c8['optional_math']['component_numerical_qualification_transferred'] is False
    assert c8['dispatch_admission']['all_observed_private_rows_have_positive_live_job']
    assert c8['dispatch_admission']['job_and_descriptor_views_match_actual_packet_arguments']
    for key, name in (('control8_acceptance', 'c8'), ('benefit', 'benefit'), ('decimal', 'decimal')):
        assert controls[key] == pin(paths[name]), key
    for label in ('a1', 'b1', 'b2', 'a2'):
        owner = json.loads((current / 'benefit/runtime' / label / 'owner.json').read_text())
        assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']


_previous_AS_authorities = authority_files
def authority_files():
    return _previous_AS_authorities() + [_AS_SOURCE / name for name in _AS_PINS]


_AT_SOURCE = Path.cwd() / 'out/native-cpu-dispatch-token-cohort/cold-frame'
_AT_PINS = {'frame_authority.py': '541b938e6f7115a5c3cf194b24795da2feb5150df190b139cf84aa16c264f26d', 'owned_guard.py': '5c92a9c3d087e16bed596d74e240d557994e87d9bd1351bc4206ecc0c757214c', 'prepare.py': '9c465d5f0b6e304d690213ec66339d15dbe7e971098a5cb85cbab3085e37d770', 'read_result.py': 'b71c168987ab5d34bf795c8f1447675ec8d359e7180573d6f83499440c804a86', 'run_owned.py': 'e2fa05f371f897ce629c60669b2fcc032fbd6b5b0ee6c5e1b63acbea87e26eb4'}
for _name, _digest in _AT_PINS.items():
    assert hashlib.sha256((_AT_SOURCE / _name).read_bytes()).hexdigest() == _digest, _name
_previous_AT_authorities = authority_files
def authority_files():
    return _previous_AT_authorities() + [_AT_SOURCE / name for name in _AT_PINS]
