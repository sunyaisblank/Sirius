"""Use read-only bounded native console geometry and derived message widths."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,ast
w=Path(__file__).resolve().parent
def load(n,h):
 b=(w/n).read_bytes();assert hashlib.sha256(b).hexdigest()==h;return b.decode()
cs=load('owned_console_running.cs','5deec7ad4156747c43a0fa4eb23642fb449aa9e64e78819110c033fd954a91a6')
cs=cs.replace('si.flags=9;si.showWindow=0;si.xCountChars=1024;si.yCountChars=2048;', 'si.flags=1;si.showWindow=0;')
cs=cs.replace(' [DllImport("kernel32.dll",SetLastError=true)] static extern bool SetConsoleScreenBufferSize(IntPtr h,Coord size);\n','')
cs=cs.replace(' public static int LastAttachError {get;private set;}', ''' public static int LastAttachError {get;private set;}
 static int minimumWidth,observedWidth,observedHeight;
 public static void ConfigureWidth(int width) {
  if(width<1 || width>4096 || launch.process!=IntPtr.Zero || attached)throw new Exception("Invalid prelaunch width contract");
  minimumWidth=width;
 }''')
cs=cs.replace('  // Geometry was requested at creation; Read verifies it without resizing history.', '  // Read-only screen access; actual geometry is bound on the first snapshot.')
old='  if(info.size.x!=1024 || info.size.y!=2048 || info.cursor.y<0 || info.cursor.y>=2047)throw new Exception("Console geometry or scroll bound invalid");'
assert cs.count(old)==1
cs=cs.replace(old,'''  if(minimumWidth<1 || info.size.x<minimumWidth || info.size.x>4096 || info.size.y<2 || info.cursor.y<0 || info.cursor.y>=info.size.y-1 || info.size.x*(info.cursor.y+1)>2097152)throw new Exception("Console geometry/read/wrap/scroll bound invalid: "+info.size.x+"x"+info.size.y);
  if(observedWidth==0){observedWidth=info.size.x;observedHeight=info.size.y;}
  if(info.size.x!=observedWidth || info.size.y!=observedHeight)throw new Exception("Console geometry changed during owned session");''')
with (w/'owned_console_native.cs').open('x') as f:f.write(cs)
ps=load('console_readiness_running.ps1','c1d17536f6b12b31a480294fdeb8e32ebe04bcc4ec60c4bc9a1b6790617ea91b')
ps=ps.replace('console-readiness-running','console-readiness-native').replace('owned_console_running.cs','owned_console_native.cs')
ps=ps.replace('[string]$LiveRevision)', '[string]$LiveRevision,[Parameter(Mandatory=$true)][int]$RequiredConsoleWidth)')
needle=' $consolePid=[SiriusOwnedConsole]::StartSuspended'
assert ps.count(needle)==1
ps=ps.replace(needle,' [SiriusOwnedConsole]::ConfigureWidth($RequiredConsoleWidth)\n $record.required_console_width=$RequiredConsoleWidth\n'+needle)
ps=ps.replace(' $record.console_geometry_requested_at_creation=@{width=1024;height=2048;startup_flags=9}\n','')
with (w/'console_readiness_native.ps1').open('x') as f:f.write(ps)
py=load('run_console_readiness_running.py','9989727dc8b7c1cf05541ee9d2b9b6259ad4922a79f0f7688aa6557e70b48e4b')
py=py.replace('console-readiness-running','console-readiness-native').replace('console_readiness_running.ps1','console_readiness_native.ps1').replace('owned_console_running.cs','owned_console_native.cs').replace('prepare_running_console.py','prepare_native_console.py')
py=py.replace("paths += [WORK/'owned_console_native.cs',WORK/'prepare_native_console.py']", "paths += [WORK/'owned_console_native.cs',WORK/'prepare_native_console.py',WORK/'capture-plan.json']")
needle="ps='/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'"
assert py.count(needle)==1
py=py.replace(needle,'''plan=read(WORK/'capture-plan.json');assert plan['source_revision']==head
science=win(ROOT/plan['stage']/'tests/backend/Release/sirius_backend_tests.exe')
future_trace=win(WORK/'native-capture-console/original-modern-te.rgp')
required_width=max(len('Connecting to process: '+science+' (PID: '+str(2**32-1)+')'),len('Trace saved to: '+future_trace),len('Profiling (RGP) feature enabled (auto-capture enabled)'),len('Waiting for GPU application to connect...'))+1
ps='/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe' ''')
py=py.replace("'-LiveRevision',head]", "'-LiveRevision',head,'-RequiredConsoleWidth',str(required_width)]")
needle="screen['width']==1024 and screen['height']==2048 and screen['cursor_y']<2047"
assert py.count(needle)==1
py=py.replace(needle,"required_width<=screen['width']<=4096 and 2<=screen['height']<=32767 and 0<=screen['cursor_y']<screen['height']-1")
needle="  assert set(screen['attached_pids'])=={o['owner_pid'],o['test_pid']}"
assert py.count(needle)==1
py=py.replace(needle,needle+"\n  initial=read(output/'initial-screen.json');assert (initial['width'],initial['height'])==(screen['width'],screen['height']) and set(initial['attached_pids'])==set(screen['attached_pids'])\n  assert o['required_console_width']==required_width and 0<len(screen['cells'])<=2097152")
ast.parse(py)
with (w/'run_console_readiness_native.py').open('x') as f:f.write(py)
for n in ('owned_console_native.cs','console_readiness_native.ps1','run_console_readiness_native.py',Path(__file__).name):
 p=w/n;print(n,p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest())
