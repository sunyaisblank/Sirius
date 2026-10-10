"""Bounded read-only terminal observer, with early Windows birth on failure."""
from pathlib import Path
import json,os,signal,subprocess,sys,time
from native_contracts import ROOT,WORK,read_case,document,sha
mode,label,case='cleanup','trial','trial'
observed={'scope':'Read-only final census of completed candidate trial, no scientific launch'}
finished=WORK
proof=ROOT/'attestations/cleanup/2026-10-10/homogeneous-opcode-dispatch-final-proofs'
assert not proof.exists();proof.mkdir(parents=True)
assert finished.is_dir()
prefix='cleanup-trial'
receipt=proof/(prefix+'.json');bridge=proof/(prefix+'-bridge.json')
temp=proof/'temp'/prefix;assert not temp.exists() and not receipt.exists() and not bridge.exists();temp.mkdir(parents=True)
def win(path):return subprocess.check_output(['wslpath','-w',str(path.resolve())],text=True).strip()
def quote(word):return "'"+word.replace("'","''")+"'"
audit=win(WORK/'trial-final-cleanup-audit.ps1');inputs=win(WORK/'trial-final-cleanup-inputs.json');bootstrap=win(temp/'bootstrap.json')
command="$ErrorActionPreference='Stop';$p=Get-Process -Id $PID;$h=$p.Handle;@{owner_pid=$PID;owner_start_time_utc=$p.StartTime.ToUniversalTime().ToString('o');expected_test_executable='';scope='Read-only terminal auditor bootstrap, no Add-Type or scientific launch'}|ConvertTo-Json|Set-Content -LiteralPath "+quote(bootstrap)+" -Encoding UTF8;& "+quote(audit)+" -Inputs "+quote(inputs)+" -Receipt "+quote(win(receipt))
arguments=['/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe','-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-Command',command]
cleanup_args=[arguments[0],'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(WORK/'emergency_native_cleanup.ps1'),'-TemporaryDirectory',win(temp),'-Output',win(temp),'-Receipt',win(proof/(prefix+'-emergency-cleanup.json'))]
record={'mode':mode,'label':label,'case':case,'scientific_read_case_invoked':False,'observer_sha256':sha(Path(__file__)),'audit_script_sha256':sha(WORK/'trial-final-cleanup-audit.ps1'),'arguments':arguments,'deadline_seconds':60,'outer_stop':False,'errors':[],'started_epoch':time.time(),'qualification_claimed':False,'scope':'Completed candidate trial saved births and provider bytes checked by bounded read-only Windows terminal audit. Hidden fields/handles/global or unconditional descendant absence not asserted.'}
started=time.monotonic()
with (proof/(prefix+'.stdout')).open('wb') as out,(proof/(prefix+'.stderr')).open('wb') as err:
 child=subprocess.Popen(arguments,stdout=out,stderr=err,start_new_session=True)
 try:
  record['linux_child_pid']=child.pid;record['linux_child_birth']=Path('/proc',str(child.pid),'stat').read_text().rsplit(')',1)[1].split()[19]
  code=child.wait(timeout=60)
 except BaseException as error:
  record['outer_stop']=True;record['errors'].append(repr(error))
 finally:
  if record['outer_stop'] or child.poll() is None:
   try:
    with (proof/(prefix+'-emergency.stdout')).open('wb') as cout,(proof/(prefix+'-emergency.stderr')).open('wb') as cerr:
     c=subprocess.run(cleanup_args,stdout=cout,stderr=cerr,timeout=60)
    assert c.returncode==0
    cleanup=document(proof/(prefix+'-emergency-cleanup.json'));assert cleanup['recorded_births_absent'] and not cleanup['errors']
   except BaseException as error:record['errors'].append(repr(error))
  try:code=child.wait(timeout=15)
  except BaseException as error:
   record['errors'].append(repr(error))
   try:os.killpg(child.pid,signal.SIGKILL)
   except ProcessLookupError:pass
   except BaseException as error:record['errors'].append(repr(error))
   try:code=child.wait(timeout=10)
   except BaseException as error:record['errors'].append(repr(error));code=child.returncode
record['returncode']=code;record['elapsed_seconds']=time.monotonic()-started
try:current=Path('/proc',str(child.pid),'stat').read_text().rsplit(')',1)[1].split()[19]
except (FileNotFoundError,ProcessLookupError):current=None
record['linux_child_birth_absent']=record.get('linux_child_birth') is not None and current!=record['linux_child_birth']
record['accepted']=code==0 and not record['outer_stop'] and not record['errors'] and record['linux_child_birth_absent']
if record['accepted']:
 terminal=document(receipt);assert terminal['pass'] and terminal['known_births_verified']==len(document(WORK/'trial-final-cleanup-inputs.json')['windows_births'])
 assert not terminal['matches'] and not terminal['accessible_compilers']
 record['terminal_receipt_sha256']=sha(receipt)
 record['bootstrap']=document(temp/'bootstrap.json')
 # A completed interop bridge normally establishes auditor exit. Keep its exact
 # birth for the next/final census, without claiming hidden consumer ownership.
bridge.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ('mode','label','case','returncode','elapsed_seconds','accepted')}),flush=True)
assert record['accepted']
result={'pass':True,'case':case,'mode':mode,'label':label,'record':observed,'terminal_receipt':str(receipt.relative_to(ROOT)),'terminal_receipt_sha256':sha(receipt),'terminal_bridge_sha256':sha(bridge),'bootstrap':record['bootstrap'],'scope':'Completed candidate trial ownership/provider cleanup acceptance with limited accessible terminal census, no dispatch or numerical/performance/full-frame/qualification claim.'}
output=proof/(prefix+'-accepted.json');assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n')
