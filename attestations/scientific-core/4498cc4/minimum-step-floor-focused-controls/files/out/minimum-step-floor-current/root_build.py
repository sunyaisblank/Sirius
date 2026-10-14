from pathlib import Path
import json,os,subprocess,sys,time
work=Path(__file__).resolve().parent
scratch=work/'tmp';scratch.mkdir()
environment=dict(os.environ,TMPDIR=str(scratch.resolve()))
with (work/'root-build-operator.log').open('xb') as log:
 for phase in ['configure','build-backend']:
  start=time.monotonic()
  result=subprocess.run(['/usr/bin/python3.12','-B',str(work/'run_owned.py'),phase],env=environment,stdout=log,stderr=subprocess.STDOUT)
  print(json.dumps({'phase':phase,'exit':result.returncode,'operator_seconds':round(time.monotonic()-start,3)}),flush=True)
  if result.returncode:sys.exit(result.returncode)
print('STRICT_CORRECTED_CONFIGURE_SOURCE_GOVERNANCE_AND_BACKEND_BUILD_COMPLETE')
