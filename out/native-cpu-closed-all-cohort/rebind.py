"""Pin exercised AW C8 authorities and bind actual AX closed3/local9 names."""
from pathlib import Path
import hashlib
_AW_SOURCE=Path.cwd()/'out/native-cpu-closed-dp-unwind-cohort'
_AW_PINS={'run_owned.py': '7d832b5c9edc5518e0aca52e1325368744ae40215bf7978f4ead4ab2740672ef', 'runtime_transition.py': 'd3c3b1a3c6115c55e791132e45c793c64c7c0c26b2aa9980382d40ae880dde19', 'optional_observer.py': '040b2656fa3e83cef5aecfb75ac12aeab6fdc8425d5e9bab4e6a0b8beaf2d68d', 'read_result.py': '94e821db4f0e14db2daf09bdc5d5a8267b3247d40c779b04a8d7f33390c63d26', 'optional_math.py': '90f322b6151d966bb53e63ca0db45e1e5c15484679bc053311165dfed997fd20', 'rebind.py': '9f5c6a09b2934d8afe093cb02faf37bb23c4f3a24b35d0f2e42d6ae87b6a8007', 'dispatch_observer.py': '1d7ac7309ff1bcf6188cb8643ef2e25de3e11dd7ecce71dacbd025747c509bc8', 'run_case.py': 'c8aa0b3f3622a278f1edc3839d4334e6c1dc61d814b68fa5425a6609dd7d20fd', 'native_observer.py': '2e6371d7a1d5fd92333ff6b9e4cb19bedd986361b3e13e519c0b6df7153b6725', 'compatibility.py': 'e983afe1c1f206bf480a41a896e64c47dd4e1e07efb0cda1b8818080b4cdedd6', 'prepare.py': '09469389c5356252fb26c1b192be4dd003bd46824e2ba4f574c1c19953a17a64', 'owned_guard.py': 'ab5a3e58e47311e363842cfa04715f475483ee05f8bd63391a4b8a67787c0d92'}
for _name,_sha in _AW_PINS.items():
    assert hashlib.sha256((_AW_SOURCE/_name).read_bytes()).hexdigest()==_sha,_name
exec(compile((_AW_SOURCE/'rebind.py').read_text(),str(_AW_SOURCE/'rebind.py'),'exec'))
_previous_AX_adapt=adapt

def adapt(text,name):
    text=_previous_AX_adapt(text,name)
    for before,after in (
        ('native-cpu-closed-dp-unwind-cohort','native-cpu-closed-all-cohort'),
        ('native-cpu-closed-dp-unwind-provider','native-cpu-closed-all-provider'),
        ('native-cpu-closed-dp-unwind-component','native-cpu-closed-all-component'),
        ('sirius_scope_transport_packet','sirius_closed_transport_packet'),
        ('sirius_scope_endpoint_packet','sirius_closed_endpoint_packet'),
        ('providerAW','providerAX'),('current_AW','current_AX')):
        text=text.replace(before,after)
    text=text.replace(
        'Only the positive-job DP v4 callback changes to a closed canonical SSA implementation with local9 cc92 strict math; unchanged AU T/E use their original scoped raw2/shared9 math.',
        'Positive genuine-job v4 callbacks select closed canonical T/E/DP implementations with one existing local9 cc92 strict math module and independently joined full unwind/GP proofs. All public/scalar/v3/AVX2 fallbacks remain exact; current original C8 supplies finite numerical evidence only within its observed paths.')
    return text

_previous_AX_authorities=authority_files
def authority_files():
    return _previous_AX_authorities()+[_AW_SOURCE/name for name in _AW_PINS]
