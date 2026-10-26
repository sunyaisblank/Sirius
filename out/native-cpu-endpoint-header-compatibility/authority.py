"""Reuse the exact AO full38 producer chain, with AP composition paths and scope."""
from pathlib import Path
import hashlib
PINS = {'prepare.py': 'd85635bdf4898ee8111572f48e025201f45c7319429e769fe463117229b11910', 'link.py': '561ce93088a2f5297b11cba8a5f2b0a64825799fe579b379a67bb48e9c2269b8', 'run_task.py': '049572331e8b5756cc9f723d7f6d4f262a0e50b339e72d40dcb737c2b3822de1', 'run_owned.py': '085bfd55e59022b1a3aae40c6f6aadcc60ebddd5a7045dd3c6e35fa04e491895', 'owned_guard.py': '246c8953a8c029c982cbcf6e8baab61487f2dd49c97d68188dee4a748465afc2'}


def load(name):
    path = Path.cwd() / 'out/native-cpu-v4-compatibility' / name
    data = path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == PINS[name], path
    text = data.decode()
    for before, after in (
        ('bin/linux-gcc/native-cpu-v4-compatibility', 'bin/linux-gcc/native-cpu-endpoint-header-compatibility'),
        ('out/native-cpu-v4-provider/provider', 'out/native-cpu-endpoint-header-provider/provider'),
        ('Reuse exact genuine R fixture driver, renamed Q baseline and unchanged AN raw3/U scalar3/pureF. Only the8th merged adapter composes exact original v3 math with separately guarded same-eight-lane optional v4 math/baseline dispatch. One new link, zero C++/raw/host/SPIRV compilations.',
         'Reuse exact genuine R fixture driver, renamed Q baseline, unchanged AN raw3 and six AO objects including the original merged guarded v4/v3 math. AP replaces only Endpoint scalar/pure admission for relative header word11241, separately checked by its focused fixture. One diagnostic link, zero diagnostic C++/raw/host/SPIRV compilations; the two component C++ objects were compiled separately.'),
        ('One fresh guarded-v4-composition Transport compatibility workload', 'One fresh AP composition Transport compatibility workload'),
        ('exact original-Q/unchanged-AN raw Transport with new guarded optional math composition', 'exact original-Q/unchanged-AN raw Transport with the AP eight-object composition')):
        text = text.replace(before, after)
    return text


def authority_files():
    root = Path.cwd()
    return ([root / 'out/native-cpu-v4-compatibility' / name for name in PINS] +
            [root / 'out/native-cpu-local-identities-compatibility' / name for name in PINS])
