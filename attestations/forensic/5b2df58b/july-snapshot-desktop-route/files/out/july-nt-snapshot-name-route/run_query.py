from pathlib import Path
import json,subprocess,sys
Q=Path(__file__).resolve().parent
O=Path(sys.argv[1]);assert O.is_dir() and not (Q/'query-bootstrap.json').exists()
PS='/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
query=[PS,'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File', '\\\\wsl.localhost\\Ubuntu\\home\\astra\\.project\\Sirius\\out\\july-nt-snapshot-name-route\\query.ps1']
audit=[PS,'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File', '\\\\wsl.localhost\\Ubuntu\\home\\astra\\.project\\Sirius\\out\\july-nt-snapshot-name-route\\terminal.ps1','-Bootstrap', '\\\\wsl.localhost\\Ubuntu\\home\\astra\\.project\\Sirius\\out\\july-nt-snapshot-name-route\\query-bootstrap.json']
record={'query_argv':query,'audit_argv':audit,'query_exit':None,'query_timeout':False,'audit_exit':None,'passed':False}
child=None
try:
 with (O/'query.stdout').open('wb') as out,(O/'query.stderr').open('wb') as err:
  child=subprocess.Popen(query,stdout=out,stderr=err,start_new_session=False)
  try:record['query_exit']=child.wait(timeout=35)
  except subprocess.TimeoutExpired:record['query_timeout']=True
 with (O/'terminal.stdout').open('wb') as out,(O/'terminal.stderr').open('wb') as err:
  terminal=subprocess.run(audit,stdout=out,stderr=err,timeout=18)
 record['audit_exit']=terminal.returncode
 if child.poll() is None:record['query_exit']=child.wait(timeout=4)
 record['passed']=record['audit_exit']==0 and not record['query_timeout'] and record['query_exit'] in (0,2)
finally:
 (O/'operations.json').write_text(json.dumps(record,indent=2)+'\n')
 if child is not None and child.poll() is None:child.terminate();child.wait(timeout=3)
assert record['passed'],record
print(json.dumps(record),flush=True)
