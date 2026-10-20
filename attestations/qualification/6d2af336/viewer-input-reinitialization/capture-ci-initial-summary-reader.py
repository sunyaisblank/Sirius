"""Capture terminal CI jobs/logs without rerunning any build or test."""
from pathlib import Path
import hashlib,json,subprocess
B=Path(__file__).resolve().parent;D=B/'ci';D.mkdir(exist_ok=True)
REPO='sunyaisblank/Sirius';REV='6d2af336b12dc430792bd5a83c04fb5632868c04'
def api(path):return subprocess.check_output(['gh','api','repos/'+REPO+'/'+path],timeout=120)
def document(path):return json.loads(api(path))
def save(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
review={};pending=False
for rid in [39487479545,39487799597]:
 run=document('actions/runs/'+str(rid));jobs=document('actions/runs/'+str(rid)+'/jobs?per_page=100')
 assert run['head_sha']==REV and run['run_attempt']==1
 assert len(jobs['jobs'])==8
 states=[]
 for job in jobs['jobs']:
  states.append({'name':job['name'],'status':job['status'],'conclusion':job['conclusion']})
  if job['status']!='completed':pending=True;continue
  if job['name'].startswith('integration-'):
   assert job['conclusion']=='success'
   path=D/(str(rid)+'-'+str(job['id'])+'.log')
   if not path.exists():path.write_bytes(api('actions/jobs/'+str(job['id'])+'/logs'))
   text=path.read_text(errors='replace')
   assert '100% tests passed, 0 tests failed out of 10' in text
   save(D/(str(rid)+'-'+str(job['id'])+'.json'),job)
  else:assert job['conclusion']=='skipped'
 if run['status']=='completed':
  assert run['conclusion']=='success'
  save(D/('run-'+str(rid)+'.json'),run);save(D/('jobs-'+str(rid)+'.json'),jobs)
  artifacts=document('actions/runs/'+str(rid)+'/artifacts?per_page=100')
  save(D/('artifacts-'+str(rid)+'.json'),artifacts)
 review[str(rid)]={'status':run['status'],'conclusion':run['conclusion'],'jobs':states}
if pending:
 print(json.dumps({'terminal':False,'runs':review}));raise SystemExit(0)
# Bind the PR's actual merge checkout separately from its source-head label.
ref=document('git/ref/pull/62/merge');commit=document('git/commits/'+ref['object']['sha']);source=document('git/commits/'+REV)
assert REV in [p['sha'] for p in commit['parents']]
assert commit['tree']['sha']==source['tree']['sha']
save(D/'PR-merge-ref.json',ref);save(D/'PR-merge-commit.json',commit);save(D/'source-commit.json',source)
files={str(p.relative_to(B)):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in D.iterdir() if p.is_file()}
result={'schema':'sirius-viewer-input-CI-review-v1','passed':True,'source_revision':REV,'runs':review,'integration_jobs_per_run':3,'original_authority_controls_per_run':30,'unrequested_full_qualification_jobs_per_run':5,'PR_merge_tree_equals_source':True,'PR_merge_commit':commit['sha'],'native_build_domains_admitted':0,'CI_artifact_bytes_downloaded_or_backed_up':False,'files':files,'scope':'Terminal exact-source and PR-merge strict complete-topology integration and original authority controls. No full Mandatory, sanitizer/scientific/provider/native-window/frame/installed-release qualification; build artifacts remain available in GitHub with API identities only.'}
save(B/'ci-review.json',result)
print(json.dumps({'terminal':True,'passed':True,'runs':list(review),'PR_merge_commit':commit['sha'],'retained_files':len(files)}))
