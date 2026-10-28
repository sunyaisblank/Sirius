"""Pin AU's original frame contract; bind only current AX evidence and scope."""
from pathlib import Path
import hashlib
import importlib.util

_AU_SOURCE = Path.cwd() / 'out/native-cpu-shared-value-cohort/cold-frame'
_AU_PINS = {
    'frame_authority.py': 'b686f15a540838a5d322723f810c7090e2667680edb4c36f018030aebe85bc2f',
    'prepare.py': '0d482987a130d4b8f9802144d6cfc5d2a2aa6da8d5150817d903b72d803d66f0',
    'read_result.py': 'b71c168987ab5d34bf795c8f1447675ec8d359e7180573d6f83499440c804a86',
    'read_terminal.py': '45333b4f90c36c9633e1fc07a848a3f15154b9f841a4b75ea6e451925819edb9',
    'run_owned.py': '64ca72c33d9057020828d8c86c75bb7983b4774175ec46dafbe5069ba594ded0',
    'owned_guard.py': '5c92a9c3d087e16bed596d74e240d557994e87d9bd1351bc4206ecc0c757214c',
}
for _name, _digest in _AU_PINS.items():
    assert hashlib.sha256((_AU_SOURCE / _name).read_bytes()).hexdigest() == _digest, _name

# These two rooted literals are the AU live-acceptance and ABBA-input roots.
# Historical AT/AU authority source paths and their byte pins stay untouched.
_source = (_AU_SOURCE / 'frame_authority.py').read_text()
assert _source.count("root / 'out/native-cpu-shared-value-cohort") == 2
_source = _source.replace("root / 'out/native-cpu-shared-value-cohort",
                          "root / 'out/native-cpu-closed-all-cohort")
exec(compile(_source, str(_AU_SOURCE / 'frame_authority.py'), 'exec'))


def _current(text):
    # Apply to generated executable source, not to historical authority pins.
    for before, after in (
        ('native-cpu-dispatch-token-cohort', 'native-cpu-closed-all-cohort'),
        ('native-cpu-shared-value-cohort', 'native-cpu-closed-all-cohort'),
        ('out/native-cpu-dispatch-token-provider/provider', 'out/native-cpu-closed-all-provider/provider'),
        ('out/native-cpu-shared-value-provider/provider', 'out/native-cpu-closed-all-provider/provider'),
        ('current_AT_link_proof', 'current_AX_link_proof'),
        ('current_AU_link_proof', 'current_AX_link_proof'),
        ('fresh AU all8/all22 controls passed', 'fresh AX all8/all22 controls passed'),
        ('AT only after all8/all22', 'AX only after all8/all22'),
        ('current AU final admission', 'current AX final admission'),
    ):
        text = text.replace(before, after)
    return text


_previous_AX_ao = ao


def ao(name):
    text = _current(_previous_AX_ao(name))
    if name == 'prepare.py':
        before = ('exact AU provider retaining AT public/scalar/pure admission and all AS public raw/GCC math/AVX2/v3 fallbacks; '
                  'only genuine positive-job v4 private callbacks enter new strict raw3 and shared-nine-core value math inside the unchanged full saved/masked/restored MXCSR scope, '
                  'same original actual-Vulkan host/scene/counts/layouts/scheduler/resources and all original numerical, control and work limits')
        after = ('exact AX provider retaining AT public/scalar/pure admission and AS public raw/GCC/AVX2/v3 fallbacks; '
                 'genuine positive-job v4 callbacks select closed canonical T/E/DP and one existing strict local-effects regcall math9 inside the complete saved/masked/restored MXCSR scope, '
                 'same original actual-Vulkan H host/scene/counts/layouts/scheduler/resources and original numerical/control/work limits; '
                 'fresh AX C8 observed E/DP SIMD/local-core PCs and Transport12 scalar, with no inferred Transport SIMD execution or prior frame acceptance')
        assert text.count(before) == 1
        text = text.replace(before, after, 1)
    return text


_previous_AX_load = load


def load(name):
    # AU's read_terminal.py repaired this exact evidence-root namespace only
    # after its failed frame. Bind the corrected current root from the outset.
    text = _current(_previous_AX_load(name))
    if name == 'read_result.py':
        assert "current = root / 'out/native-cpu-closed-all-cohort'" in text
    return text


_previous_AX_authorities = authority_files


def authority_files():
    return _previous_AX_authorities() + [_AU_SOURCE / name for name in _AU_PINS]


def _current_c8():
    path = Path.cwd() / 'out/native-cpu-closed-all-cohort/benefit/correctness.py'
    assert pin(path) == dict(bytes=8266, sha256='9cd6469007d85f73971dc1944f73e3f2339e5d6f6fc60a7dbe6ffd4dbd0b7712')
    spec = importlib.util.spec_from_file_location('AX_frame_current_C8_gate', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return set(module.check()) | {path}


_previous_AX_acceptance_inputs = acceptance_inputs


def acceptance_inputs():
    return _previous_AX_acceptance_inputs() | _current_c8()


_previous_AX_check_acceptance = check_acceptance


def check_acceptance():
    """Require current AX all8/all22/Decimal/C8/static evidence before any mode."""
    _previous_AX_check_acceptance()
    _current_c8()
    root = Path.cwd()
    current = root / 'out/native-cpu-closed-all-cohort'
    controls = json.loads((current / 'controls/readers/controls-1-8-acceptance.json').read_text())
    c8 = json.loads((current / 'control-8-acceptance.json').read_text())
    assert controls['reader'] == pin(current / 'controls/readers/read_controls.py')
    assert controls['optional_C8_runtime_coverage'] == c8['optional_math']
    assert set(controls['remaining7']) == {'control-' + str(i) for i in range(1, 8)}
    actions = json.loads((current / 'controls/runtime/actions.json').read_text())
    for label, record in controls['remaining7'].items():
        folder = current / 'controls/runtime' / label
        owner = json.loads((folder / 'owner.json').read_text())
        assert record['passed'] and record['owner'] == pin(folder / 'owner.json')
        assert record['XML'] == pin(folder / 'gtest.xml')
        assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
        assert owner['source'] == dict(revision=actions[label]['revision'], status='')
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
        assert owner['observed_births'] and all(item['absent'] for item in owner['observed_births'])
        assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
        assert owner['limits'] == actions[label]['limits'] == dict(seconds=300, rss_kib=4194304, reserve_kib=2097152)
