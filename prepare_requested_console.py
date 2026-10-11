"""Request only the derived required width at creation, with read-only verification."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,ast
w=Path(__file__).resolve().parent
def load(n,h):
 b=(w/n).read_bytes();assert hashlib.sha256(b).hexdigest()==h;return b.decode()
cs=load('owned_console_native.cs','f527cd3ab4f1ae1f2a4be85414240108797ec1d0018664161c11fb058457dac0')
needle='si.flags=1;si.showWindow=0;'
assert cs.count(needle)==1
cs=cs.replace(needle,'si.flags=9;si.showWindow=0;si.xCountChars=(uint)minimumWidth;si.yCountChars=2048;')
with (w/'owned_console_requested.cs').open('x') as f:f.write(cs)
ps=load('console_readiness_native.ps1','7af03a3825817b0f7056f2e10b432df552febac600b3f886b0433039d31a4ed8')
ps=ps.replace('console-readiness-native','console-readiness-requested').replace('owned_console_native.cs','owned_console_requested.cs')
with (w/'console_readiness_requested.ps1').open('x') as f:f.write(ps)
py=load('run_console_readiness_native.py','876dedf300cafc5061e2a90286d275e35541c8c3a4fa137f7ebe819b16c9008a')
py=py.replace('console-readiness-native','console-readiness-requested').replace('console_readiness_native.ps1','console_readiness_requested.ps1').replace('owned_console_native.cs','owned_console_requested.cs').replace('prepare_native_console.py','prepare_requested_console.py')
ast.parse(py)
with (w/'run_console_readiness_requested.py').open('x') as f:f.write(py)
for n in ('owned_console_requested.cs','console_readiness_requested.ps1','run_console_readiness_requested.py',Path(__file__).name):
 p=w/n;print(n,p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest())
