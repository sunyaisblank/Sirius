"""Original finite 4B observer, with one engine-policy observation addon."""
from pathlib import Path
import hashlib

P = Path(__file__).resolve().parents[2] / 'lavapipe-acqrel-ordering/ordinary/capture.py'
assert hashlib.sha256(P.read_bytes()).hexdigest() == '803758e9ce81b34330283c9e81a99bff92647ad176e60a10e1679086b82a40dd'
source = P.read_text()
changes = {
    'bin/linux-gcc/lavapipe-acqrel-ordering/src/gallium/targets/lavapipe/libvulkan_lvp.so':
        'bin/linux-gcc/native-cpu-arithmetic-sharing/provider/libvulkan_lvp.so',
    'out/native-original-input-ownership/compiler-result.json':
        'out/native-cpu-arithmetic-sharing/compiler-result.json',
    "env=os.environ.copy();": "exec(compile((W/'symbol_bindings.py').read_text(), str(W/'symbol_bindings.py'), 'exec'))\nenv=os.environ.copy();",
    "env.update(SIRIUS_GEOMETRY_CAPTURE=": "env.update(**codegen_environment);env.update(SIRIUS_GEOMETRY_CAPTURE=",
    "str(W/'properties_capture.py')+'\\nrun\\n'":
        "str(W/'properties_capture.py')+'\\nsource '+str(W/'capture_codegen.py')+'\\nrun\\n'",
    'Exact prior admitted complete module, new guarded acquire-release provider':
        'Exact unchanged cooperative T and newly shared native CPU E/DP, scoped complete bulk-zero Default provider;4Bguard only',
}
for old, new in changes.items():
    assert source.count(old) == 1, old
    source = source.replace(old, new)
# __file__ intentionally remains this wrapper, so all observer/output paths
# remain in this namespace. The underlying immutable probe/preload are reused.
exec(compile(source, str(P), 'exec'))
