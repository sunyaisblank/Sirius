"""Observe actual initial console geometry, then size only the owned console."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,ast
w=Path(__file__).resolve().parent
def load(n,h):
 b=(w/n).read_bytes();assert hashlib.sha256(b).hexdigest()==h;return b.decode()
cs=load('owned_console_running.cs','5deec7ad4156747c43a0fa4eb23642fb449aa9e64e78819110c033fd954a91a6')
old=' public static Snapshot Read(uint owner,uint child) {\n'
assert cs.count(old)==1
cs=cs.replace(old,''' public static Snapshot Read(uint owner,uint child){return Capture(owner,child,true);}
 public static Snapshot Observe(uint owner,uint child){return Capture(owner,child,false);}
 public static void SizeOwnedScreen() {
  if(!attached || screen==IntPtr.Zero)throw new Exception("No owned console attachment");
  if(!SetConsoleScreenBufferSize(screen,new Coord(1024,2048)))throw Error();
 }
 static Snapshot Capture(uint owner,uint child,bool strict) {
''')
old='  if(info.size.x!=1024 || info.size.y!=2048 || info.cursor.y<0 || info.cursor.y>=2047)throw new Exception("Console geometry or scroll bound invalid");'
assert cs.count(old)==1
cs=cs.replace(old,'''  if(info.size.x<1 || info.size.x>4096 || info.size.y<1 || info.cursor.y<0 || info.cursor.y>=info.size.y || info.size.x*(info.cursor.y+1)>2097152)throw new Exception("Initial console observation bound invalid");
  if(strict && (info.size.x!=1024 || info.size.y!=2048 || info.cursor.y>=2047))throw new Exception("Console geometry or scroll bound invalid: "+info.size.x+"x"+info.size.y+" cursor "+info.cursor.x+","+info.cursor.y);''')
cs=cs.replace("if(line[line.Length-1]!=' ')","if(strict && line[line.Length-1]!=' ')")
with (w/'owned_console_sized.cs').open('x') as f:f.write(cs)
ps=load('console_readiness_running.ps1','c1d17536f6b12b31a480294fdeb8e32ebe04bcc4ec60c4bc9a1b6790617ea91b')
ps=ps.replace('console-readiness-running','console-readiness-sized').replace('owned_console_running.cs','owned_console_sized.cs')
old=' $initial=[SiriusOwnedConsole]::Read([uint32]$PID,[uint32]$process.Id)'
assert ps.count(old)==1
ps=ps.replace(old,''' $before=[SiriusOwnedConsole]::Observe([uint32]$PID,[uint32]$process.Id)
 $before | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'pre-sizing-screen.json') -Encoding UTF8
 [SiriusOwnedConsole]::SizeOwnedScreen()
 $record.console_sized_before_any_scientific_application=$true
 $initial=[SiriusOwnedConsole]::Read([uint32]$PID,[uint32]$process.Id)''')
with (w/'console_readiness_sized.ps1').open('x') as f:f.write(ps)
py=load('run_console_readiness_running.py','9989727dc8b7c1cf05541ee9d2b9b6259ad4922a79f0f7688aa6557e70b48e4b')
py=py.replace('console-readiness-running','console-readiness-sized').replace('console_readiness_running.ps1','console_readiness_sized.ps1').replace('owned_console_running.cs','owned_console_sized.cs').replace('prepare_running_console.py','prepare_sized_console.py')
py=py.replace("o['console_geometry_verified_after_attachment']", "o['console_geometry_verified_after_attachment'] and o['console_sized_before_any_scientific_application']")
needle="  assert set(screen['attached_pids'])=={o['owner_pid'],o['test_pid']}"
assert py.count(needle)==1
py=py.replace(needle,needle+"\n  before=read(output/'pre-sizing-screen.json');assert set(before['attached_pids'])=={o['owner_pid'],o['test_pid']} and 0<len(before['cells'])<=2097152\n  assert len(before['cells'])==before['width']*before['rows_read']")
ast.parse(py)
with (w/'run_console_readiness_sized.py').open('x') as f:f.write(py)
for n in ('owned_console_sized.cs','console_readiness_sized.ps1','run_console_readiness_sized.py',Path(__file__).name):
 p=w/n;print(n,p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest())
