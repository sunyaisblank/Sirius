"""Pin the exercised AO wrappers and original AN authorities; bind fresh AP paths."""
from pathlib import Path
import hashlib
AO_PINS = {'prepare.py': 'a97f3d86d1a61f31390215a2a8b84981ae451337408277f8cdd86e077f7d003f', 'run_case.py': '987242b691936c347986d2c6de2cfb014c8545e4c84b39280714d136c1020f83', 'read_result.py': '23c8a22bf8bc0faf5133d479a1fe48ef68f5b7d2480fc45b2940a50f98b01dde', 'native_observer.py': 'a99b8d4dda4c5f21fc711a66af7772dfbc92e884916b0a5d1f30792ef646917f', 'compatibility.py': '551446b585794ee7d5a3f608b9551cd07ee0afc37d1b0692226abf36cfe034be', 'runtime_transition.py': 'e21d338c1b7173a193dfd845e2daacfe4e42c2ced058c341940e5f13f31863ea', 'optional_math.py': '749d23a7fe11e8e1f738e524f6fa7de3203b9858ad29d83d834acf997f64eace', 'optional_observer.py': '3bf880a5d1be6d9710f03dec6637a992f7e985f083681e66c2846711cab65334', 'run_owned.py': 'f2da572ee2f392bb5b076c86216aa05a07ff7fb61cd6c5e179b5a2a0ab8d5c11', 'owned_guard.py': 'a8a9bdafd2340c56372032b549c3e3766cadb4dcc8bdb36657e877348e089a5b'}
PINS = {'prepare.py': '8b7c0c099cafc68ade3a6f1fbdb61efc3f1e505a7e2e3f2767da1a567584bb2d', 'run_case.py': '70a89563785e82d28a43574a8914740df2b6045c90d19e360b0fd9b1714642cf', 'read_result.py': '5903eff3ad76cd652b7447ca76cf112131bcc0717608282729e804fcb0a60837', 'native_observer.py': 'd1e55e7eed62836c425e369ee43e865d69b71530e6f046adba80f33526ca085b', 'runtime_transition.py': '1afe9523fa43639ad5fc9e14609c98f024a1a8368bccec3dc10a92e1bf150f3e', 'run_owned.py': '085bfd55e59022b1a3aae40c6f6aadcc60ebddd5a7045dd3c6e35fa04e491895', 'owned_guard.py': '246c8953a8c029c982cbcf6e8baab61487f2dd49c97d68188dee4a748465afc2'}


def ao(name):
    path = Path.cwd() / 'out/native-cpu-v4-cohort' / name
    data = path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == AO_PINS[name], path
    text = data.decode()
    for before, after in (
        ('out/native-cpu-v4-cohort', 'out/native-cpu-endpoint-header-cohort'),
        ('out/native-cpu-v4-provider/provider', 'out/native-cpu-endpoint-header-provider/provider'),
        ('out/native-cpu-v4-compatibility', 'out/native-cpu-endpoint-header-compatibility'),
        ('providerAO', 'providerAP'), ('current_AO', 'current_AP'),
        ('Fresh AO guarded-math combined Transport compatibility', 'Fresh AP composition Transport compatibility'),
        ('current AO provider object-code association', 'current AP provider object-code association'),
        ('Same original H host/AO guarded arithmetic provider', 'Same original H host/AP Endpoint admission provider')):
        text = text.replace(before, after)
    before = 'AO changes only private object composition/identity, preserving actual existing v3 eligibility, FP/precision properties, raw3/scheduler/fallback and original physical contracts.'
    after = 'AP changes only Endpoint pure/scalar admission for relative header word11241 and provider identity; the unchanged shader retains its acceptance/refusal and all output behavior. Other six AO objects, actual v3 eligibility, FP/precision properties, raw3/scheduler/fallback and original physical contracts remain.'
    return text.replace(before, after)


def load(name):
    path = Path.cwd() / 'out/native-cpu-local-identities-cohort' / name
    data = path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == PINS[name], path
    text = data.decode().replace('out/native-cpu-local-identities-cohort', 'out/native-cpu-endpoint-header-cohort')
    text = text.replace('out/native-cpu-local-identities/provider', 'out/native-cpu-endpoint-header-provider/provider')
    text = text.replace('providerAN', 'providerAP').replace('current_AN', 'current_AP')
    return text.replace('identical U scalar3', 'unchanged Transport/DopriPhase scalar objects and verified Endpoint scalar/pure admission extension at relative word11241')


def authority_files():
    root = Path.cwd()
    return ([root / 'out/native-cpu-v4-cohort' / name for name in AO_PINS] +
            [root / 'out/native-cpu-local-identities-cohort' / name for name in PINS])


def seal_authorities(text):
    anchor = 'files = set(expected)'
    assert text.count(anchor) == 1
    return text.replace(anchor, 'for path in authority_files():\n    expected[str(path)] = pin(path)\n' + anchor, 1)
