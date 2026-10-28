"""Pin the exercised AP producer chain; bind only AS sources and raw symbols."""
from pathlib import Path
import hashlib

_AR_SOURCE = Path.cwd() / 'out/native-cpu-packed-math-cohort'
_AR_PINS = {'compatibility.py': 'f2ae0fb85f8e6646fdee2a271bb30b0eb9b0ca38391ce4ca263989f6fd96b488', 'native_observer.py': 'd351a37aded9549738e05548ed1e3a983d3d50de6256b324dc0a7f298d559190', 'optional_math.py': 'd3fd4d943fffd2f83534c7e3582bf045ead2c86aa6e44bc6b9f4f140669d041b', 'optional_observer.py': '6526a8455c67b3330093129394d8cf740c5caaeaaf7502f05b8afb4514e17523', 'owned_guard.py': 'ab5a3e58e47311e363842cfa04715f475483ee05f8bd63391a4b8a67787c0d92', 'prepare.py': '118cce14e369f5844c9435d5e85d5a63e67926df6ec38028a96c62884060d28d', 'read_result.py': '94e821db4f0e14db2daf09bdc5d5a8267b3247d40c779b04a8d7f33390c63d26', 'rebind.py': 'c45fed40589d386f11721b931e155e291a6dae6567a6a758a3a174dda6dd78aa', 'run_case.py': 'c8aa0b3f3622a278f1edc3839d4334e6c1dc61d814b68fa5425a6609dd7d20fd', 'run_owned.py': '7d832b5c9edc5518e0aca52e1325368744ae40215bf7978f4ead4ab2740672ef', 'runtime_transition.py': 'd3c3b1a3c6115c55e791132e45c793c64c7c0c26b2aa9980382d40ae880dde19'}
for _name, _digest in _AR_PINS.items():
    assert hashlib.sha256((_AR_SOURCE / _name).read_bytes()).hexdigest() == _digest, _name

R = Path.cwd()
AP = R / 'out/native-cpu-endpoint-header-cohort'
AP_PINS = {'run_owned.py': '7d832b5c9edc5518e0aca52e1325368744ae40215bf7978f4ead4ab2740672ef', 'runtime_transition.py': 'd3c3b1a3c6115c55e791132e45c793c64c7c0c26b2aa9980382d40ae880dde19', 'optional_observer.py': '6526a8455c67b3330093129394d8cf740c5caaeaaf7502f05b8afb4514e17523', 'read_result.py': '94e821db4f0e14db2daf09bdc5d5a8267b3247d40c779b04a8d7f33390c63d26', 'optional_math.py': '873b7eed86d2234210d146d80369ac9cd622dc031b545dbc6daaf34f7f84ca5f', 'rebind.py': '420331da24d7e5768532d2466a6a1d7b80036579204c38aa2f82ef0af4371c95', 'run_case.py': 'c8aa0b3f3622a278f1edc3839d4334e6c1dc61d814b68fa5425a6609dd7d20fd', 'native_observer.py': 'e3819bb87b6133f30d98cb1a4bc4623a7bedf5d9e41de8d715383722c191fccf', 'compatibility.py': '2d8b51030d408bca810ce94b586f07258113fca12b45cfc0a9872a0116b0c224', 'prepare.py': '118cce14e369f5844c9435d5e85d5a63e67926df6ec38028a96c62884060d28d', 'owned_guard.py': 'ab5a3e58e47311e363842cfa04715f475483ee05f8bd63391a4b8a67787c0d92'}
for name, digest in AP_PINS.items():
    assert hashlib.sha256((AP / name).read_bytes()).hexdigest() == digest, name
_AP = {}
exec(compile((AP / 'rebind.py').read_text(), str(AP / 'rebind.py'), 'exec'), _AP)


def adapt(text, name):
    for before, after in (
        ('out/native-cpu-endpoint-header-cohort', 'out/native-cpu-inline-private-cohort'),
        ('out/native-cpu-endpoint-header-provider/provider', 'out/native-cpu-inline-private-provider/provider'),
        ('out/native-cpu-endpoint-header-compatibility', 'out/native-cpu-v4-raw-packets'),
        ('providerAP', 'providerAS'), ('current_AP', 'current_AS'),
        ('Same original H host/AP Endpoint admission provider', 'Same original H host/AS guarded raw provider')):
        text = text.replace(before, after)
    # Preserve actual observed function names; adjust only the new symbol prefix
    # and inverse stage/parent mapping, never normalize a fallback observation.
    text = text.replace("'sirius_simd_' + slug + '_packet'", "'sirius_simd_v4_' + slug + '_packet'")
    text = text.replace("len('sirius_simd_')", "len('sirius_simd_v4_')")
    text = text.replace(".replace('sirius_simd_', 'sirius_native_')", ".replace('sirius_simd_v4_', 'sirius_native_')")
    if name == 'optional_observer.py':
        before = "for index, slug in enumerate(('transport', 'endpoint', 'dopri_phase'), 4):\n    item = proof['objects'][index]"
        after = "for slug, item in zip(('transport', 'endpoint', 'dopri_phase'), proof['raw_v4_objects']):"
        assert text.count(before) == 1
        text = text.replace(before, after, 1)
    if name == 'optional_math.py':
        before = "raw = dict(zip(('transport', 'endpoint', 'dopri_phase'), proof['objects'][4:7]))"
        after = "raw = dict(zip(('transport', 'endpoint', 'dopri_phase'), proof['raw_v4_objects']))"
        assert text.count(before) == 1
        text = text.replace(before, after, 1)
    text = text.replace(
        'AP changes only Endpoint pure/scalar admission for relative header word11241 and provider identity; the unchanged shader retains its acceptance/refusal and all output behavior. Other six AO objects, actual v3 eligibility, FP/precision properties, raw3/scheduler/fallback and original physical contracts remain.',
        'AS preserves AP first7 native objects/admission/scalar and original AVX2 fallback. Its eighth object combines same-eight-lane raw-v4 stages with forced opcode inlining, zeroed privateD and576B compiler-private math banks, using the unchanged reviewed zero-residual arithmetic behind the unchanged baseline OS-aware packet selector. All9380 operators, original FP/layout/counts/scheduler/fallback/physical contracts remain; new raw-v4 paths require current runtime observation.')
    return text


def ao(name):
    return adapt(_AP['ao'](name), name)


def load(name):
    return adapt(_AP['load'](name), name)


def authority_files():
    return _AP['authority_files']() + [AP / name for name in AP_PINS] + [_AR_SOURCE / name for name in _AR_PINS]


def seal_authorities(text):
    anchor = 'files = set(expected)'
    assert text.count(anchor) == 1
    return text.replace(anchor, 'for path in authority_files():\n    expected[str(path)] = pin(path)\n' + anchor, 1)
