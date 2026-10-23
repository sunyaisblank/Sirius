"""Prepare one original four-byte fallback capture for each native stage."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
O = W / 'ordinary'
O.mkdir(exist_ok=False)


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


reuse = json.loads((R / 'out/native-cpu-arithmetic-sharing/ordinary/reuse.json').read_text())
copied = {}
for name in ('capture_jit.py', 'capture_features.py', 'properties_capture.py',
             'read_features.py', 'float-layout.json', 'feature-layout.json'):
    item = reuse['observer_sources'][name]
    source = R / item['resolved_regular_source']
    assert pin(source) == {key: item[key] for key in ('bytes', 'sha256')}
    (O / name).write_bytes(source.read_bytes())
    copied[name] = {'source': str(source.relative_to(R)), **pin(source)}
(O / 'native_observer.py').write_bytes((W / 'native_observer.py').read_bytes())
base = R / reuse['core_capture']['path']
assert pin(base)['sha256'] == reuse['core_capture']['sha256']
wrapper = '''from pathlib import Path
import hashlib
source_path = Path.cwd() / 'BASE'
assert hashlib.sha256(source_path.read_bytes()).hexdigest() == 'DIGEST'
source = source_path.read_text()
changes = CHANGES
for old, new in changes.items():
    assert source.count(old) == 1, old
    source = source.replace(old, new)
exec(compile(source, str(source_path), 'exec'))
'''
changes = {
    'bin/linux-gcc/lavapipe-acqrel-ordering/src/gallium/targets/lavapipe/libvulkan_lvp.so':
        'bin/linux-gcc/native-cpu-slang-bridge/provider/libvulkan_lvp.so',
    'out/native-original-input-ownership/compiler-result.json':
        'out/native-cpu-row-ownership/compiler-result.json',
    "str(W/'properties_capture.py')+'\\nrun\\n'":
        "str(W/'properties_capture.py')+'\\nsource '+str(W/'native_observer.py')+'\\nrun\\n'",
    "print(json.dumps({'passed':True,'module_bytes':loaded['bytes'],'total_JIT_function_bytes':sum(x['bytes'] for x in event['functions'])}))":
        "native=json.loads((folder/'native-calls.json').read_text()); assert len(native)==1 and native[0]['returned'] and not native[0]['handled']; assert native[0]['output_before']==native[0]['output_after']; print(json.dumps({'passed':True,'module_bytes':loaded['bytes'],'native_fallback_before_writes':True}))",
    'Exact prior admitted complete module, new guarded acquire-release provider':
        'Exact original row module and newly linked native CPU bridge; original4B negative guard falls back before writes',
}
for old in changes:
    assert base.read_text().count(old) == 1, old
(O / 'capture.py').write_text(wrapper.replace('BASE', str(base.relative_to(R)))
                             .replace('DIGEST', pin(base)['sha256']).replace('CHANGES', repr(changes)))
for name in ('run_owned.py', 'owned_guard.py'):
    (O / name).write_bytes((W / 'build' / name).read_bytes())
(O / 'reuse.json').write_text(json.dumps({'core': reuse['core_capture'], 'observers': copied,
    'scope': 'Original4B guard/probe and passive observers, exact new provider + original all-row module paths; add passive native entry/return/output snapshot. No inferior calls or changed workloads.'}, indent=2) + '\n')
print(json.dumps({'prepared_guard_stages': 3, 'executed': False}))
