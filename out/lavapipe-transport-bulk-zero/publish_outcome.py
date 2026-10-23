"""Publish the exact preserved rejection and retain complete readback receipts."""
from pathlib import Path
import hashlib,json,subprocess,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent
D=R/'attestations/qualification/06adfcdc/lavapipe-transport-bulk-zero-finite-outcome'
COMMIT=subprocess.check_output(['git','rev-parse','8e028435'],text=True).strip()
def doc(p):return json.loads(p.read_text())
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
def api(path,*args):return json.loads(subprocess.check_output(['gh','api',path,*args],text=True))
m=doc(D/'MANIFEST.json');f=doc(D/'findings.json')
assert f['benefit_checks_passed']==18 and not f['all22_original_benefit_gates_pass']
assert not f['all8_original_controls_admitted'] and not f['cold_frame']['attempted']
assert all(f['original_controls'][str(i)]=={'status':'not_run'} for i in range(1,8))
for name,value in m['archives'].items():assert pin(D/name)==value
remote=subprocess.check_output(['git','ls-remote','origin','refs/heads/development/lavapipe-transport-bulk-zero-2026-10-23'],text=True).split()[0]
assert subprocess.run(['git','merge-base','--is-ancestor',COMMIT,remote]).returncode==0
assert subprocess.check_output(['git','show',COMMIT+':'+str((D/'findings.json').relative_to(R))])==(D/'findings.json').read_bytes()
link='https://github.com/sunyaisblank/Sirius/tree/'+COMMIT+'/'+str(D.relative_to(R))
previous='https://github.com/sunyaisblank/Sirius/issues/106#issuecomment-6356912234'
paragraph=f'''The latest private complete-row Transport bulk-zero extension at source06adfcdc is **rejected:18/22** (2/4whole-call,14/16inclusive-stage,2/2provisioning-inclusiveAPI; independent Decimal agrees). Original C8 passed20.892055s with all32complete readbacks, and all four sole-ABBA slots passed correctness/audit. Forward/reverse legacy whole-call gains were -5.230190%/-3.210851% versus required+5%; legacy Transport gains -30.640668%/-25.311536% failed the original -5% floor. Modern whole-call gains18.968556%/27.653368% and both API gains above85% do not override those failures. Remaining7 and the frame were not run; no retry or adoption. [Preserved source, raw evidence and actual products]({link}).

This distinct extension applied the exact complete73440-byte clear to the immutable full row-T hash/type for the first time, retaining all fields/inputs/masks/lifetime and original Default optimization. One incremental provider build3.098441s reused858sealed inputs; one6.166356s4B guard capture observed the resolved-libc clear and allzero return, with actual Default/greedy and full remaining-IR equality under193bijectively renamed SSA locals. No Boolean/arithmetic rewrite or capture retry. Current T evidence composes disclosed historical unchanged E/DP source/code; no fresh E/DP capture claim. Authentic7f2 all-row host/archive/three SPVs were restored with zero new Sirius/shader/probe compilation. Seven owned runs are terminal/reaped. A dual-T host route was rejected source-only because it would break the original7kernel-load assertion and lacked a proved crossover; no cohort was spent on it.

The [preceding mixed cooperative-T/bulk-E/DP candidate]({previous}) still supplies22/22 and all8 finite results, but its single original frame stopped300.235003s at2827/16384 without completedXML. Both outcomes remain private finite evidence. The original frame, full science, July image identity, Mandatory, native-window/Radeon scenes and exact-main installed-release readiness remain open.

'''
saved=doc(W/'publication-issue44-before.json');current=api('repos/sunyaisblank/Sirius/issues/44')
assert current['body']==saved['body'],'Issue changed; reconcile before writing.'
start=current['body'].index('The latest private bulk-zero candidate at source337d0bae')
end_marker='Issue44 and the ultimate goal remain open.\n\n'
end=current['body'].index(end_marker,start)+len(end_marker)
body=current['body'][:start]+paragraph+current['body'][end:]
assert len(body)<65000
payload=W/'publication-issue44-payload.json';payload.write_text(json.dumps({'body':body})+'\n')
updated=api('repos/sunyaisblank/Sirius/issues/44','--method','PATCH','--input',str(payload));assert updated['body']==body
readback=api('repos/sunyaisblank/Sirius/issues/44');assert readback['body']==body
(W/'publication-issue44-readback.json').write_text(json.dumps(readback,indent=2)+'\n')
comment='Complete Transport bulk-zero finite outcome: '+paragraph+'\nAll original gates and public source0109ab46 remain unchanged; issue44 and the ultimate goal remain open.'
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
