"""One original finite guard capture, selecting the immutable full Transport row module."""
from pathlib import Path
import hashlib
_basis=Path.cwd()/'out/lavapipe-owner-bulk-zero/ordinary/capture.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest()=='3d703bac7d812c90fa5dbe301ca8fa929111d2ed6b30587c91ac628a95ce1e12'
_source=_basis.read_text()
for _old,_new in {
 'bin/linux-gcc/lavapipe-owner-bulk-zero/provider/libvulkan_lvp.so':'bin/linux-gcc/lavapipe-transport-bulk-zero/provider/libvulkan_lvp.so',
 'out/native-cooperative-transport-intergroup/compiler-result.json':'out/native-cpu-row-ownership/compiler-result.json',
 'Exact mixed cooperative Transport/row Endpoint/Dopri modules, scoped complete bulk-zero Default provider':'Exact full row Transport module, scoped complete bulk-zero Default provider',
}.items():
 assert _source.count(_old)==1;_source=_source.replace(_old,_new)
exec(compile(_source,str(_basis),'exec'))
