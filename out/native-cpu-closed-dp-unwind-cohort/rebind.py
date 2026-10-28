"""Reuse every pinned AV C8 assertion/observer; bind only the AW repair."""
from pathlib import Path
import hashlib
_AV_SOURCE=Path.cwd()/'out/native-cpu-closed-dp-cohort'
_AV_PINS={'run_owned.py': '7d832b5c9edc5518e0aca52e1325368744ae40215bf7978f4ead4ab2740672ef', 'runtime_transition.py': 'd3c3b1a3c6115c55e791132e45c793c64c7c0c26b2aa9980382d40ae880dde19', 'optional_observer.py': 'eb323f83e4952271a3dfd3c97473e07bdbd6dff3b7c701f6660850dbe831d2f3', 'read_result.py': '94e821db4f0e14db2daf09bdc5d5a8267b3247d40c779b04a8d7f33390c63d26', 'optional_math.py': '9b04321c73fc1eb275ee14237aa8385fd97cdc026c1368013e6f0b0fbb2975bf', 'rebind.py': '54dca4e075e4b26800b7a6769cfb3942a116e518d4a438a5631b56baac9b9ccd', 'dispatch_observer.py': '4d41f408c9b15eb1f5017dd3f655e9ed65b36e657d3f5ac6991c28c355a78124', 'run_case.py': 'c8aa0b3f3622a278f1edc3839d4334e6c1dc61d814b68fa5425a6609dd7d20fd', 'native_observer.py': 'ed886853b1f4bb3cd16374fdbaedc4e1b3fb856a9a58f3f1fe9419c649442b37', 'compatibility.py': '28f4357f697f2fe60a0f4e17e8fb1cdd04e022feccc8f52ec6bc87e04835850b', 'prepare.py': 'fab1d8ed8358faa24d492424aca7e5a7f9d9feb648fee871fdafc412412ad830', 'owned_guard.py': 'ab5a3e58e47311e363842cfa04715f475483ee05f8bd63391a4b8a67787c0d92'}
for _name,_sha in _AV_PINS.items():
    assert hashlib.sha256((_AV_SOURCE/_name).read_bytes()).hexdigest()==_sha,_name
exec(compile((_AV_SOURCE/'rebind.py').read_text(),str(_AV_SOURCE/'rebind.py'),'exec'))
_previous_AW_adapt=adapt

def adapt(text,name):
    text=_previous_AW_adapt(text,name)
    for before,after in (
        ('native-cpu-closed-dp-cohort','native-cpu-closed-dp-unwind-cohort'),
        ('native-cpu-closed-dp-provider','native-cpu-closed-dp-unwind-provider'),
        ('native-cpu-closed-dp-component','native-cpu-closed-dp-unwind-component'),
        ('providerAV','providerAW'),('current_AV','current_AW')):
        text=text.replace(before,after)
    return text

_previous_AW_authorities=authority_files
def authority_files():
    return _previous_AW_authorities()+[_AV_SOURCE/name for name in _AV_PINS]
