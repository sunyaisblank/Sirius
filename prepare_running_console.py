"""Correct suspended-console readiness using creation-time console dimensions."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,ast
w=Path(__file__).resolve().parent
def load(n,h):
 b=(w/n).read_bytes();assert hashlib.sha256(b).hexdigest()==h;return b.decode()
cs=load('owned_console.cs','ea9e3dc969b685fc3939fc0b8cdba7d178990b29209e383399fa00fdad7d27f9')
cs=cs.replace('si.flags=1;si.showWindow=0;', 'si.flags=9;si.showWindow=0;si.xCountChars=1024;si.yCountChars=2048;')
cs=cs.replace('ResumeThread(launch.thread)==UInt32.MaxValue','ResumeThread(launch.thread)!=1')
cs=cs.replace(' static bool attached=false;', ' static bool attached=false;\n public static int LastAttachError {get;private set;}')
cs=cs.replace('if(!AttachConsole(pid))return false;', 'if(!AttachConsole(pid)){LastAttachError=Marshal.GetLastWin32Error();return false;}\n  LastAttachError=0;')
cs=cs.replace('  // Wide enough for bounded identity/save lines; bottom-row use is rejected.\n  if(!SetConsoleScreenBufferSize(screen,new Coord(1024,2048)))throw Error();', '  // Geometry was requested at creation; Read verifies it without resizing history.')
with (w/'owned_console_running.cs').open('x') as f:f.write(cs)
ps=load('console_readiness.ps1','d3ac3fb3da8cfc7408b3ecdec2811a182b4f112c1cc60f081cff5f5aace32317')
ps=ps.replace('console-readiness','console-readiness-running').replace('owned_console.cs','owned_console_running.cs')
old=" $attached=$false;$record.console_ready_observed=$false\n $attach="
assert ps.count(old)==1
ps=ps.replace(old," $record.console_geometry_requested_at_creation=@{width=1024;height=2048;startup_flags=9}\n [SiriusOwnedConsole]::Resume()\n $attached=$false;$record.console_ready_observed=$false\n $attach=")
ps=ps.replace('Owned suspended-console attachment/lifetime/RSS bound exceeded','Owned running-console attachment/lifetime/RSS bound exceeded')
needle='  $attached=[SiriusOwnedConsole]::TryAttach([uint32]$process.Id)\n  if (!$attached)'
assert ps.count(needle)==1
ps=ps.replace(needle,'  $attached=[SiriusOwnedConsole]::TryAttach([uint32]$process.Id)\n  $record.last_console_attach_error=[SiriusOwnedConsole]::LastAttachError\n  if (!$attached)')
ps=ps.replace('$record.console_geometry_initialized_before_resume=$true','$record.console_geometry_verified_after_attachment=$true')
needle=" $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8\n [SiriusOwnedConsole]::Resume()\n while (!$process.HasExited)"
assert ps.count(needle)==1
ps=ps.replace(needle," $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8\n while (!$process.HasExited)")
with (w/'console_readiness_running.ps1').open('x') as f:f.write(ps)
py=load('run_console_readiness.py','beee10a456440529c1fdc62d05c062c8a09e70c4b653102eec66a8d668027775')
py=py.replace('console-readiness','console-readiness-running').replace('console_readiness.ps1','console_readiness_running.ps1').replace('owned_console.cs','owned_console_running.cs').replace('prepare_console_readiness.py','prepare_running_console.py').replace('console_geometry_initialized_before_resume','console_geometry_verified_after_attachment')
ast.parse(py)
with (w/'run_console_readiness_running.py').open('x') as f:f.write(py)
for n in ('owned_console_running.cs','console_readiness_running.ps1','run_console_readiness_running.py',Path(__file__).name):
 p=w/n;print(n,p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest())
