"""Reuse sealed AN C8 authorities with current artifact paths only."""
from pathlib import Path
import hashlib
PINS = {'prepare.py': '8b7c0c099cafc68ade3a6f1fbdb61efc3f1e505a7e2e3f2767da1a567584bb2d', 'run_case.py': '70a89563785e82d28a43574a8914740df2b6045c90d19e360b0fd9b1714642cf', 'read_result.py': '5903eff3ad76cd652b7447ca76cf112131bcc0717608282729e804fcb0a60837', 'native_observer.py': 'd1e55e7eed62836c425e369ee43e865d69b71530e6f046adba80f33526ca085b', 'runtime_transition.py': '1afe9523fa43639ad5fc9e14609c98f024a1a8368bccec3dc10a92e1bf150f3e', 'run_owned.py': '085bfd55e59022b1a3aae40c6f6aadcc60ebddd5a7045dd3c6e35fa04e491895', 'owned_guard.py': '246c8953a8c029c982cbcf6e8baab61487f2dd49c97d68188dee4a748465afc2'}
def load(name):
    p=Path.cwd()/'out/native-cpu-local-identities-cohort'/name
    data=p.read_bytes()
    assert hashlib.sha256(data).hexdigest()==PINS[name],name
    source=data.decode().replace('out/native-cpu-local-identities-cohort','out/native-cpu-v4-cohort')
    source=source.replace('out/native-cpu-local-identities/provider','out/native-cpu-v4-provider/provider')
    source=source.replace('providerAN','providerAO').replace('current_AN','current_AO')
    return source
