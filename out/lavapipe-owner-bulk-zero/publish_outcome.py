"""Publish only preserved finite results, with exact remote readback receipts."""
from pathlib import Path
import hashlib,json,subprocess,sys,zipfile
R=Path.cwd(); W=Path(__file__).resolve().parent
D=R/'attestations/qualification/337d0bae/lavapipe-owner-bulk-zero-finite-outcome'
COMMIT=subprocess.check_output(['git','rev-parse',sys.argv[1]],text=True).strip()
BRANCH='development/lavapipe-owner-bulk-zero-2026-10-23'
def doc(p):return json.loads(p.read_text())
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
def api(path,*args):return json.loads(subprocess.check_output(['gh','api',path,*args],text=True))
manifest=doc(D/'MANIFEST.json');findings=doc(D/'findings.json');frame=findings['cold_frame']
for name,value in manifest['archives'].items():assert pin(D/name)==value
remote=subprocess.check_output(['git','ls-remote','origin','refs/heads/'+BRANCH],text=True).split()[0]
assert subprocess.run(['git','merge-base','--is-ancestor',COMMIT,remote]).returncode==0
assert subprocess.check_output(['git','show',COMMIT+':'+str((D/'findings.json').relative_to(R))])==(D/'findings.json').read_bytes()
assert findings['all22_original_benefit_gates_pass'] and findings['all8_original_controls_admitted']
assert findings['independent_Decimal_agrees'] and frame['attempted']
link='https://github.com/sunyaisblank/Sirius/tree/'+COMMIT+'/'+str(D.relative_to(R))
gains='/'.join(f"{r['improvement_percent']:.6f}%" for r in findings['whole_call_gains_percent'])
detail=frame['disposition'];progress=detail.get('last_progress')
frame_text=(f"The single original cold Kerr frame {'passed' if frame['accepted'] else 'did not pass'}: "
 f"{frame['wall_seconds']:.6f}s, stop={frame['stop_reason']}, "
 f"last observed progress={progress}, completed XML={frame['XML_present']}.")
paragraph=f'''The latest private bulk-zero candidate at source337d0bae passed **all22 original benefit checks and all8 original controls**. Independent Decimal agrees; original C8 passed21.165241s with all32complete readbacks, and both precision routes passed28complete readbacks in control1. Sole A1/B1/B2/A2,3internal repeats/slot: whole-call gains forward legacy/modern and reverse legacy/modern were {gains}; all16inclusive stages met the original -5% floor and both provisioning-inclusive API gains exceeded85%. No pooling, subtraction, resampling, warmup, retry or threshold change.

{frame_text} The original scene, cache-disable policy, 300s/4GiB/2GiB owner, render assertions and governors were retained. Process ownership has been closed; owner success alone is not frame acceptance. [Preserved source, raw evidence and actual products]({link}).

The scoped change keeps original Default machine optimization and replaces only exact admitted Endpoint/DopriPhase array initialization with an equivalent complete bulk zero at the same point/lifetime. Actual E/DP emitted code and resolved libc binding observed every byte zero; remaining optimized IR equivalence includes explicitly reviewed defined-mask Boolean rewrites, not a machine-checked proof. Cooperative Transport is excluded from this source gate and reuses its disclosed historical source/code proof; no new Tcapture is claimed. One incremental provider build21.857724s reused739sealed inputs; the mixed f958 host and all3SPVs were reused with zero new Sirius/shader/probe compilation, followed by one original-object frame link. These finite results do not admit full science, original July image identity, native-window/Radeon scenes, Mandatory, exact-main installation or release readiness. Issue44 and the ultimate goal remain open.

'''
current=api('repos/sunyaisblank/Sirius/issues/44')
(W/'publication-issue44-current-before.json').write_text(json.dumps(current,indent=2)+'\n')
first,rest=current['body'].split('\n\n',1);body=first+'\n\n'+paragraph+rest
need=api('repos/sunyaisblank/Sirius/issues/44');assert need['body']==current['body'],'Issue changed; reconcile before writing.'
payload=W/'publication-issue44-payload.json';payload.write_text(json.dumps({'body':body})+'\n')
updated=api('repos/sunyaisblank/Sirius/issues/44','--method','PATCH','--input',str(payload));assert updated['body']==body
readback=api('repos/sunyaisblank/Sirius/issues/44');assert readback['body']==body
(W/'publication-issue44-readback.json').write_text(json.dumps(readback,indent=2)+'\n')
comment='Bulk-zero finite outcome: '+paragraph
payload=W/'publication-issue106-payload.json';payload.write_text(json.dumps({'body':comment})+'\n')
created=api('repos/sunyaisblank/Sirius/issues/106/comments','--input',str(payload));assert created['body']==comment
receipt=api('repos/sunyaisblank/Sirius/issues/comments/'+str(created['id']));assert receipt['body']==comment
(W/'publication-issue106-readback.json').write_text(json.dumps(receipt,indent=2)+'\n')
files=[W/'publish_outcome.py',*[p for p in W.glob('publication-*.json')],W/'preregistration.txt',W/'preregistration-payload.json',W/'preregistration-readback.json']
archive=D/'publication-receipts.zip';assert not archive.exists()
with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,str(p.relative_to(R)))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for p in files:assert z.read(str(p.relative_to(R)))==p.read_bytes()
record={'source_preserved_at_commit':COMMIT,'issue44_exact_body_readback':True,'issue106_exact_comment_readback':True,
 'issue106_comment':receipt['html_url'],'fixed_protocol_preregistration':doc(W/'preregistration-readback.json')['html_url'],
 'receipt_archive':pin(archive),'receipt_members':{str(p.relative_to(R)):pin(p) for p in files},'full_goal_complete':False}
(D/'publication.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'issue44_updated':True,'issue106':receipt['html_url'],'publication_receipts':pin(archive)}))
