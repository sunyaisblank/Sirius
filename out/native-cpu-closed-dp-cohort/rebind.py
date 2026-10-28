"""Reuse the pinned AU original C8 authorities; bind the mixed AV implementation."""
from pathlib import Path
import hashlib
_AU_SOURCE=Path.cwd()/'out/native-cpu-shared-value-cohort'
_AU_PINS={'run_owned.py': '7d832b5c9edc5518e0aca52e1325368744ae40215bf7978f4ead4ab2740672ef', 'runtime_transition.py': 'd3c3b1a3c6115c55e791132e45c793c64c7c0c26b2aa9980382d40ae880dde19', 'optional_observer.py': '79df71b9ca439b869f85294e27c872a34bb90b4fc649d2121a4f319d73b516ac', 'read_result.py': '94e821db4f0e14db2daf09bdc5d5a8267b3247d40c779b04a8d7f33390c63d26', 'optional_math.py': '7a40bdca1a9282edca0397b65880db6bf258659d4dcb40e1288fbd393c342aa5', 'rebind.py': 'ad37e3a271da3f74a5fbcc98fad09e61b062a1243be2a3094a4776c4b31a95ca', 'dispatch_observer.py': 'cba16eb61b2086ffac92587c64eb53ceaf23170b8a3407804b37baa90ceb8cf9', 'run_case.py': 'c8aa0b3f3622a278f1edc3839d4334e6c1dc61d814b68fa5425a6609dd7d20fd', 'native_observer.py': '43078bed96efa7a41a93aec9c871a6f626d45867021c495a8aa9e3628c9a055c', 'compatibility.py': '3ade549d77fb4125d556b7cdee0b0f6d6b374ef243827ea389ca9d8f2c62f6f3', 'prepare.py': 'd657613b09ad6d37a199a00edd8491ef19fbf90d19c1901698c0ca523f4fccee', 'owned_guard.py': 'ab5a3e58e47311e363842cfa04715f475483ee05f8bd63391a4b8a67787c0d92'}
for _name,_sha in _AU_PINS.items():
    assert hashlib.sha256((_AU_SOURCE/_name).read_bytes()).hexdigest()==_sha,_name
exec(compile((_AU_SOURCE/'rebind.py').read_text(),str(_AU_SOURCE/'rebind.py'),'exec'))
_previous_AV_adapt=adapt


def adapt(text,name):
    text=_previous_AV_adapt(text,name)
    for before,after in (
        ('native-cpu-shared-value-cohort','native-cpu-closed-dp-cohort'),
        ('native-cpu-shared-value-provider','native-cpu-closed-dp-provider'),
        ('native-cpu-shared-value-component','native-cpu-closed-dp-component'),
        ('sirius_scope_dp_packet','sirius_closed_dp_packet'),
        ('providerAU','providerAV'),('current_AU','current_AV'),
        ("('sirius_simd_', 'sirius_scope_')","('sirius_simd_', 'sirius_scope_', 'sirius_closed_')")):
        text=text.replace(before,after)
    text=text.replace('Only positive-job admitted v4 branches enter separately named strict raw3 and one shared nine-core value module.',
        'Only the positive-job DP v4 callback changes to a closed canonical SSA implementation with local9 cc92 strict math; unchanged AU T/E use their original scoped raw2/shared9 math.')
    return text


_previous_AV_authorities=authority_files
def authority_files():
    return _previous_AV_authorities()+[_AU_SOURCE/name for name in _AU_PINS]
