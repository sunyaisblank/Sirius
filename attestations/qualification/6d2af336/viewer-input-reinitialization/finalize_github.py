from pathlib import Path
import json,subprocess
R=Path.cwd();B=Path(__file__).resolve().parent
assert json.loads((B/'ci-review.json').read_text())['passed']
head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
assert subprocess.check_output(['git','ls-remote','--heads','origin','development/viewer-input-evidence-2026-10-21'],text=True).split()==[head,'refs/heads/development/viewer-input-evidence-2026-10-21']
link='https://github.com/sunyaisblank/Sirius/tree/'+head+'/attestations/qualification/6d2af336/viewer-input-reinitialization'
needle='[Current platform CI](https://github.com/sunyaisblank/Sirius/actions/runs/39487479545) is still running.'
replacement='[Current platform CI](https://github.com/sunyaisblank/Sirius/actions/runs/39487479545) passes all three strict integration jobs and their 30 original authority controls. [Preserved counterexample, current local products, source/build/input controls and CI records]('+link+'). Current local alignment still admits 0/8 external domains; full Mandatory/scientific/native-window/frame/installed-release qualification remains open.'
records={}
for issue in [44,48,53]:
 value=json.loads(subprocess.check_output(['gh','issue','view',str(issue),'--json','body,state,url']))
 assert value['state']=='OPEN' and needle in value['body']
 body=value['body'].replace(needle,replacement,1)
 path=B/('github-issue-'+str(issue)+'-terminal.md');path.write_text(body)
 subprocess.run(['gh','issue','edit',str(issue),'--body-file',str(path)],check=True,capture_output=True)
 current=json.loads(subprocess.check_output(['gh','issue','view',str(issue),'--json','body,state,url']))
 assert current['state']=='OPEN' and current['body']==body;records[str(issue)]=current
value=json.loads(subprocess.check_output(['gh','pr','view','62','--json','body,headRefOid,isDraft,state,url']))
assert value['headRefOid']=='6d2af336b12dc430792bd5a83c04fb5632868c04' and value['isDraft'] and value['state']=='OPEN' and needle in value['body']
body=value['body'].replace(needle,replacement,1);p=B/'github-PR62-terminal.md';p.write_text(body)
subprocess.run(['gh','pr','edit','62','--body-file',str(p)],check=True,capture_output=True)
current=json.loads(subprocess.check_output(['gh','pr','view','62','--json','body,headRefOid,isDraft,state,url']))
assert current['body']==body and current['isDraft'] and current['state']=='OPEN';records['PR62']=current
value=json.loads(subprocess.check_output(['gh','issue','view','103','--json','body,url']))
terminal='''Completed at clean pushed source [6d2af336](https://github.com/sunyaisblank/Sirius/commit/6d2af336b12dc430792bd5a83c04fb5632868c04). Assigning the complete validated CameraState clears all six held movement flags on successful reinitialisation. The original input control now checks each held key across Stop/Initialise, unchanged idle coordinates/no restart, and fresh input recovery; both strict GCC14 and Clang21 runs pass one case with zero failures/errors/skips/disabled.\n\nBoth existing strict C++26 presets build all seven roles and source/alignment gates. All 35 retained SPVs and three generated math headers remain byte-exact. CI39487479545 passes all three platform integration jobs/30 authority controls; the draft PR's CI39487799597 also passes, with actual merge tree equal to source6d2af336. Five full-qualification/diagnostic jobs per run were unrequested.\n\n'''+('[Complete preserved source/counterexample/products/owner/CI records]('+link+'). This closes the input-reset defect only. No native input-delivery, frame-speed, full Mandatory/scientific/provider/release qualification follows; current external-domain admission remains0/8. PR62 stays draft, main stays5e706073 and #44/#48/#53 remain open.\n\nOriginal finding and acceptance retained:\n\n')
p=B/'github-issue-103-terminal.md';p.write_text(terminal+value['body'])
subprocess.run(['gh','issue','edit','103','--body-file',str(p)],check=True,capture_output=True)
subprocess.run(['gh','issue','close','103','--reason','completed'],check=True,capture_output=True)
current=json.loads(subprocess.check_output(['gh','issue','view','103','--json','body,state,stateReason,url']))
assert current['state']=='CLOSED' and current['body']==terminal+value['body'];records['103']=current
(B/'github-terminal-review.json').write_text(json.dumps({'passed':True,'source_revision':'6d2af336b12dc430792bd5a83c04fb5632868c04','immutable_evidence_link':link,'records':records,'original_acceptance_preserved':True},indent=2)+'\n')
print(json.dumps({'passed':True,'issue103':'closed completed','PR62':'draft open','issues44_48_53':'open','source':'6d2af336'}))
