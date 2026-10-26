"""Retain original AN control authorities; bind current cohort paths only."""
from pathlib import Path
import hashlib
PINS={'prepare.py': '34b2ba450669d79adef3f4addd038b7e445782800b889855620efec7adc05085', 'run_controls.py': 'a347e191ee83311c843f0021ad2b9de41bd1d6f87e41ab0e150211ec6f72883c', 'readers/read_controls.py': '56308b92841eb5cb17626cfcab65a74d16ce49922257000378be8c2888c445ab', 'readers/native_evidence.py': 'b9024736eaffc834151c9b99e4e3a9b82664da274ef1b98ede7d0221e7692857', 'runtime/run_case.py': 'af700fb2c5394bb079756ae9bfa7371f26010ce2349e9d1a1161e598d8c4a354', 'runtime/run_owned.py': '6ae4510d2a1672ffb6b1fee3e37bdbc62c3aee3e993cafb1c210224b55cf2dd2', 'runtime/owned_guard.py': '6fc32fd7fbed8e24370e929829a70b82aa53e54668fb38ded68ace2d1e9a6037'}
def load(name):
    path=Path.cwd()/'out/native-cpu-local-identities-cohort/controls'/name
    data=path.read_bytes()
    assert hashlib.sha256(data).hexdigest()==PINS[name],name
    return data.decode().replace('out/native-cpu-local-identities-cohort','out/native-cpu-endpoint-header-cohort').replace('H host/AN','H host/AP').replace('host/AN','host/AP')

AO_PINS = {'prepare.py': '0dc7811dca48c90a4f62718b75f1bac231512cd3c5c590777ba8f2374ddaff0a', 'run_controls.py': '310f168a7433e224ca992db7056389fa6b2c4079c4d77727e6c90f52a26c65f0', 'readers/read_controls.py': '75c30fb2877475f7f7a926ece90216de038ed3559263a91df2506cb2d350adf1', 'readers/native_evidence.py': 'b9c382e16c8563e20f162f1523235377f07e7d787e2124a7d1115ccd9f7bdc1e', 'runtime/run_case.py': '736fc5723750cc5cd2254417effe71c9dc36bc81948a0a5eac6bf44aaed27d65', 'runtime/run_owned.py': 'eaaf90cbd4efd5790921bc2c039e03da83f759f2db60e8f72e0b1cc3414d2280', 'runtime/owned_guard.py': 'b15ba893c5fdab7b5c6abfec3bf239ec7e145f3ada66ab8e98c4cfdf158136e7'}


def ao(name):
    path = Path.cwd() / 'out/native-cpu-v4-cohort/controls' / name
    data = path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == AO_PINS[name], path
    return data.decode().replace('out/native-cpu-v4-cohort', 'out/native-cpu-endpoint-header-cohort').replace('H host/AO', 'H host/AP').replace('host/AO', 'host/AP').replace('guarded arithmetic provider', 'Endpoint header admission provider')


def authority_files():
    root = Path.cwd()
    return ([root / 'out/native-cpu-v4-cohort/controls' / name for name in AO_PINS] +
            [root / 'out/native-cpu-local-identities-cohort/controls' / name for name in PINS])
