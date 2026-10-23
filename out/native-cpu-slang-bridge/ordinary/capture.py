from pathlib import Path
import hashlib
source_path = Path.cwd() / 'out/lavapipe-acqrel-ordering/ordinary/capture.py'
assert hashlib.sha256(source_path.read_bytes()).hexdigest() == '803758e9ce81b34330283c9e81a99bff92647ad176e60a10e1679086b82a40dd'
source = source_path.read_text()
changes = {'bin/linux-gcc/lavapipe-acqrel-ordering/src/gallium/targets/lavapipe/libvulkan_lvp.so': 'bin/linux-gcc/native-cpu-slang-bridge/provider/libvulkan_lvp.so', 'out/native-original-input-ownership/compiler-result.json': 'out/native-cpu-row-ownership/compiler-result.json', "str(W/'properties_capture.py')+'\\nrun\\n'": "str(W/'properties_capture.py')+'\\nsource '+str(W/'native_observer.py')+'\\nrun\\n'", "print(json.dumps({'passed':True,'module_bytes':loaded['bytes'],'total_JIT_function_bytes':sum(x['bytes'] for x in event['functions'])}))": "native=json.loads((folder/'native-calls.json').read_text()); assert len(native)==1 and native[0]['returned'] and not native[0]['handled']; assert native[0]['output_before']==native[0]['output_after']; print(json.dumps({'passed':True,'module_bytes':loaded['bytes'],'native_fallback_before_writes':True}))", 'Exact prior admitted complete module, new guarded acquire-release provider': 'Exact original row module and newly linked native CPU bridge; original4B negative guard falls back before writes'}
for old, new in changes.items():
    assert source.count(old) == 1, old
    source = source.replace(old, new)
exec(compile(source, str(source_path), 'exec'))
