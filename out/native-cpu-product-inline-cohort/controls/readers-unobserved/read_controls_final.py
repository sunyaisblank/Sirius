"""Reconcile two nonexecuted authoring files against exact historical Git bytes."""
from pathlib import Path
import hashlib
import json
import subprocess

R = Path.cwd()
reviewed_nonexecuted_sources = {}
for name in ('out/native-cpu-product-inline-cohort/cold-frame/prepare.py',
             'out/preserve-native-cpu-product-inline.py'):
    pins = {}
    for label, revision in (('before', '5c3b3ea771694865ad4ba74237e7ac69b50b7ea4'),
                            ('after', '8499357d2c3724e0eb0e83b5dabf02cb5deb2ee7')):
        data = subprocess.check_output(['git', 'show', revision + ':' + name])
        pins[label] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
    data = (R / name).read_bytes()
    assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == pins['after']
    reviewed_nonexecuted_sources[name] = pins

basis = Path(__file__).resolve().with_name('read_controls.py')
wrapper = basis.read_text()
marker = "exec(compile(source, str(basis), 'exec'))"
assert wrapper.count(marker) == 1
injection = '''before = "        assert pin(p) == {k: value[k] for k in ['bytes', 'sha256']}, key"
after = """        expected_pin = {k: value[k] for k in ['bytes', 'sha256']}
        current_pin = pin(p)
        if current_pin != expected_pin:
            relative = str(p.relative_to(R)) if p.is_relative_to(R) else str(p)
            assert i <= 3 and relative in reviewed_nonexecuted_sources, key
            transition = reviewed_nonexecuted_sources[relative]
            assert expected_pin == transition['before'] and current_pin == transition['after']
            historical = subprocess.check_output(['git', 'show', owner['source']['revision'] + ':' + relative])
            assert {'bytes': len(historical), 'sha256': hashlib.sha256(historical).hexdigest()} == expected_pin
        else:
            assert current_pin == expected_pin, key"""
assert source.count(before) == 1
source = source.replace(before, after)
source = source.replace("'remaining7': accepted,", "'reviewed_nonexecuted_source_transitions': reviewed_nonexecuted_sources, 'remaining7': accepted,")
'''
wrapper = wrapper.replace(marker, injection + marker)
exec(compile(wrapper, str(basis), 'exec'))
