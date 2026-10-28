"""Pin AU's original seven/all8 authorities; rebind only current AX paths."""
from pathlib import Path
import hashlib

_AU_SOURCE = Path.cwd() / 'out/native-cpu-shared-value-cohort/controls'
_AU_PINS = {
    'control_authority.py': '1bbdc1b46115c8547db7edd0ef6495507204f092f91fb5ebe2f30c8ab3324d37',
    'prepare.py': 'ea86922097d90fda3100ae8acd4df925a665f19be26ccd7b4da6cfba71ae4b87',
    'run_controls.py': '97732438c0acb7b32d78e1e860d8101ac59ff8a9b9e0884e5c1e66f4f7b948cf',
    'readers/read_controls.py': 'c4bf7c87de1f3594fae8562bff830586779434e93f7e213165bb0338d98290b2',
    'readers/native_evidence.py': '45f9d7034c9fa38d25be86d0ea21d445b9f5b0cf25c91985d5858436cca4898c',
    'runtime/run_case.py': '60c7e027ca92667c339151672ee8d2982b09125af81656d997ad2bdf30a92cb0',
    'runtime/run_owned.py': '0d2204e3850f96f2aa86e6f8bb43bccd2c64a724cf7645d5e9f75fd79f8f7bfc',
    'runtime/owned_guard.py': '36573bfbb46b8f41d4fcdfdcecfb8c6475423931f7e1fc1cbd5dfe1701b6411e',
}
for _name, _digest in _AU_PINS.items():
    assert hashlib.sha256((_AU_SOURCE / _name).read_bytes()).hexdigest() == _digest, _name
exec(compile((_AU_SOURCE / 'control_authority.py').read_text(),
             str(_AU_SOURCE / 'control_authority.py'), 'exec'))

_previous_AX_load = load
_previous_AX_ao = ao


def _current(text):
    return (text.replace('out/native-cpu-shared-value-cohort', 'out/native-cpu-closed-all-cohort')
            .replace('out/native-cpu-shared-value-provider/provider', 'out/native-cpu-closed-all-provider/provider')
            .replace('H host/AU', 'H host/AX').replace('host/AU', 'host/AX'))


def load(name):
    return _current(_previous_AX_load(name))


def ao(name):
    return _current(_previous_AX_ao(name))


_previous_AX_authorities = authority_files


def authority_files():
    return _previous_AX_authorities() + [_AU_SOURCE / name for name in _AU_PINS]
