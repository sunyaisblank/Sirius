"""Retain original AN control authorities; bind current cohort paths only."""
from pathlib import Path
import hashlib
PINS={'prepare.py': '34b2ba450669d79adef3f4addd038b7e445782800b889855620efec7adc05085', 'run_controls.py': 'a347e191ee83311c843f0021ad2b9de41bd1d6f87e41ab0e150211ec6f72883c', 'readers/read_controls.py': '56308b92841eb5cb17626cfcab65a74d16ce49922257000378be8c2888c445ab', 'readers/native_evidence.py': 'b9024736eaffc834151c9b99e4e3a9b82664da274ef1b98ede7d0221e7692857', 'runtime/run_case.py': 'af700fb2c5394bb079756ae9bfa7371f26010ce2349e9d1a1161e598d8c4a354', 'runtime/run_owned.py': '6ae4510d2a1672ffb6b1fee3e37bdbc62c3aee3e993cafb1c210224b55cf2dd2', 'runtime/owned_guard.py': '6fc32fd7fbed8e24370e929829a70b82aa53e54668fb38ded68ace2d1e9a6037'}
def load(name):
    path=Path.cwd()/'out/native-cpu-local-identities-cohort/controls'/name
    data=path.read_bytes()
    assert hashlib.sha256(data).hexdigest()==PINS[name],name
    return data.decode().replace('out/native-cpu-local-identities-cohort','out/native-cpu-v4-cohort').replace('H host/AN','H host/AO').replace('host/AN','host/AO')
