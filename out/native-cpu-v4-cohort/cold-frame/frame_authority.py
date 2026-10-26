"""Reuse exercised AN cold-frame source with current AO identities only."""
from pathlib import Path
import hashlib
PINS={'prepare.py': '40a2ce8a89fde25eed6e931357ceeeb0f8d754af5c3aef0bc58ab3f23a2025f8', 'read_result.py': '0ddfc56360c5395981a8b9e8a14aadd6bd241e26dc28f7da778e3a02cefa80e9', 'run_owned.py': '085bfd55e59022b1a3aae40c6f6aadcc60ebddd5a7045dd3c6e35fa04e491895', 'owned_guard.py': '246c8953a8c029c982cbcf6e8baab61487f2dd49c97d68188dee4a748465afc2'}
def load(name):
    path=Path.cwd()/'out/native-cpu-local-identities-cohort/cold-frame'/name
    data=path.read_bytes()
    assert hashlib.sha256(data).hexdigest()==PINS[name],name
    source=data.decode().replace('native-cpu-local-identities-cohort','native-cpu-v4-cohort')
    source=source.replace('out/native-cpu-local-identities/provider','out/native-cpu-v4-provider/provider')
    source=source.replace('current_AN_link_proof','current_AO_link_proof')
    source=source.replace('fresh AN all8/all22 controls passed','fresh AO all8/all22 controls passed')
    return source
