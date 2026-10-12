"""Launch one exact non-promoting diagnostic CI export and record the request."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import datetime,hashlib,json,re,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
source=json.loads((WORK/'source.json').read_text());head=source['source_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert subprocess.check_output(['git','rev-parse','@{u}'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'],text=True)
assert json.loads((WORK/'host-source-review.json').read_text())['accepted']
request=WORK/'workflow-dispatch-request.json';response=WORK/'workflow-dispatch-response.json'
assert not request.exists() and not response.exists()
argv=['gh','workflow','run','ci.yml','--repo','sunyaisblank/Sirius','--ref',source['branch'],'-f','full_qualification=false']
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
request.write_text(json.dumps({'source_revision':head,'arguments':argv,'created_utc':now(),'scope':'Exact diagnostic source/build/non-render CI and fresh Windows native producer export only; full qualification false. No original scientific/runtime/frame/release gate weakened.'},indent=2)+'\n')
r=subprocess.run(argv,capture_output=True,text=True,timeout=60)
response.write_text(json.dumps({'source_revision':head,'arguments':argv,'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr,'response_utc':now()},indent=2)+'\n')
assert r.returncode==0
urls=re.findall(r'https://github.com/sunyaisblank/Sirius/actions/runs/(\d+)',r.stdout)
assert len(urls)==1, 'Inspect saved dispatch response before any further launch'
run={'source_revision':head,'branch':source['branch'],'run_id':int(urls[0]),'request_sha256':hashlib.sha256(request.read_bytes()).hexdigest()}
(WORK/'ci-run.json').write_text(json.dumps(run,indent=2)+'\n')
print(json.dumps(run))
