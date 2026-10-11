"""Parse the capture controller without running it; retain owned parser births."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, os, subprocess, time
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def win(p):return subprocess.check_output(['wslpath','-w',str(p.resolve())],text=True).strip()
def quote(s):return "'"+s.replace("'","''")+"'"
def birth(pid):
 try:return Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[19]
 except (FileNotFoundError,ProcessLookupError):return None
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
ps='/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
output=WORK/'counter-controls-syntax';assert not output.exists();output.mkdir();temp=output/'temp'
paths=[Path(__file__),WORK/'profiler_counter_controls.ps1',WORK/'emergency_native_cleanup.ps1']
inputs={str(p.relative_to(ROOT)):seal(p) for p in paths}
target=WORK/'profiler_counter_controls.ps1'
command="$ErrorActionPreference='Stop';$t=(New-Item -ItemType Directory -Path "+quote(win(temp))+" -Force).FullName;$p=Get-Process -Id $PID;$h=$p.Handle;@{owner_pid=$PID;owner_start_time_utc=$p.StartTime.ToUniversalTime().ToString('o');expected_test_executable='';scope='Native Parser.ParseFile only; no controller execution or child compiler'}|ConvertTo-Json|Set-Content -LiteralPath (Join-Path $t 'bootstrap.json') -Encoding UTF8;$tokens=$null;$errors=$null;$ast=[System.Management.Automation.Language.Parser]::ParseFile("+quote(win(target))+",[ref]$tokens,[ref]$errors);@{source_sha256=(Get-FileHash -LiteralPath "+quote(win(target))+" -Algorithm SHA256).Hash.ToLowerInvariant();tokens=$tokens.Count;errors=@($errors|ForEach-Object {$_.Message});passed=($errors.Count -eq 0)}|ConvertTo-Json -Depth 6|Set-Content -LiteralPath "+quote(win(output/'result.json'))+" -Encoding UTF8;if($errors.Count){exit 1}"
args=[ps,'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-Command',command]
record={'input_seals':inputs,'linux_owner_pid':os.getpid(),'linux_owner_start_ticks':birth(os.getpid()),'scope':'Native syntax parsing and exact-birth observer closure only; no device/GUI/compiler/scientific launch.','errors':[]}
def run(name,args,timeout):
 with (output/(name+'.stdout')).open('xb') as out,(output/(name+'.stderr')).open('xb') as err:
  p=subprocess.Popen(args,stdout=out,stderr=err,start_new_session=True);ticks=birth(p.pid)
  r={'arguments':args,'pid':p.pid,'start_ticks':ticks,'timeout_seconds':timeout};record[name]=r
  try:r['returncode']=p.wait(timeout=timeout)
  except BaseException:
   r['timeout_or_interruption']=True
   raise
  finally:r['linux_birth_absent']=birth(p.pid)!=ticks
 return r['returncode']
start=time.monotonic();code=None
try:
 code=run('parser',args,20)
 assert code==0 and read(output/'result.json')['passed']
 assert read(output/'result.json')['source_sha256']==inputs[str(target.relative_to(ROOT))]['sha256']
except BaseException as e:record['errors'].append(repr(e))
finally:
 # The reviewed helper targets only the parser's exact bootstrap birth and scoped temporary compiler fields.
 cleanup=[ps,'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(WORK/'emergency_native_cleanup.ps1'),'-TemporaryDirectory',win(temp),'-Output',win(output),'-Receipt',win(output/'closure.json')]
 try:
  assert run('closure',cleanup,60)==0
  closure=read(output/'closure.json');assert closure['recorded_births_absent'] and not closure['errors']
  if code==0:assert not any(t['same_birth_live'] for t in closure['targets'])
 except BaseException as e:record['errors'].append(repr(e))
record['elapsed_seconds']=time.monotonic()-start
record['accepted']=code==0 and not record['errors'] and all(v.get('linux_birth_absent') for k,v in record.items() if k in ('parser','closure'))
assert all(seal(ROOT/p)==v for p,v in inputs.items())
(output/'bridge.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'accepted':record['accepted'],'elapsed_seconds':record['elapsed_seconds'],'errors':record['errors']}),flush=True)
assert record['accepted']
