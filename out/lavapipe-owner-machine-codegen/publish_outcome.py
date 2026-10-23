"""Publish the preserved finite rejection with exact readback receipts."""
from pathlib import Path
import hashlib, json, subprocess, zipfile
R=Path.cwd(); W=Path(__file__).resolve().parent
D=R/'attestations/qualification/33623035/lavapipe-owner-machine-codegen-finite-outcome'
COMMIT='1a9d3004'
def doc(p):return json.loads(p.read_text())
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
def api(path,*args):return json.loads(subprocess.check_output(['gh','api',path,*args],text=True))
manifest=doc(D/'MANIFEST.json')
for name,value in manifest['archives'].items():assert pin(D/name)==value
remote=subprocess.check_output(['git','ls-remote','origin','refs/heads/development/lavapipe-owner-machine-codegen-2026-10-23'],text=True).split()[0]
assert subprocess.run(['git','merge-base','--is-ancestor',COMMIT,remote]).returncode==0
saved=doc(W/'publication-issue44-before.json'); current=api('repos/sunyaisblank/Sirius/issues/44')
assert current['body']==saved['body'],'Issue changed; reconcile before writing.'
link='https://github.com/sunyaisblank/Sirius/tree/'+subprocess.check_output(['git','rev-parse',COMMIT],text=True).strip()+'/attestations/qualification/33623035/lavapipe-owner-machine-codegen-finite-outcome'
paragraph=f'''The latest private scoped machine-codegen trial at source33623035 is **rejected:6/22** (0/4whole-call,4/16inclusive-stage,2/2provisioning-inclusiveAPI; independent Decimal agrees). Original C8 passed18.465985s and all32full readbacks matched; all four sole-ABBA slots also passed correctness/audit. Whole-call gains were -23.189192%/-92.653648% forward legacy/modern and -18.644068%/-78.795587% reverse. Remaining7 controls and the original frame were not run. No unchanged retry or public adoption. [Preserved source, raw evidence and actual products]({link}).

This distinct mechanism selected per-engine machine CodeGenNone only after original optimized LLVM IR for exact admitted row Endpoint/DopriPhase; actual raw optimized IR matched prior bytes, both used the fast allocator, and cooperative Transport retained greedy/no-bridge compilation. The original mixed f958 ELF and all3SPVs were reused, with one incremental provider build21.573927s and739sealed reusable build inputs unchanged. Actual E/DP code grew to2102153/442539bytes. Shorter complete-process wall time did not satisfy retained-call/stage performance. Ten owned runs are terminal/reaped; observer defects and the source-composed Transport Default premise are disclosed. This neither establishes the prior300s timeout cause nor admits science, Mandatory, frame, native-window or installed-release readiness.

'''
first,rest=current['body'].split('\n\n',1); body=first+'\n\n'+paragraph+rest
payload=W/'publication-issue44-payload.json';payload.write_text(json.dumps({'body':body})+'\n')
updated=api('repos/sunyaisblank/Sirius/issues/44','--method','PATCH','--input',str(payload));assert updated['body']==body
readback=api('repos/sunyaisblank/Sirius/issues/44');assert readback['body']==body
(W/'publication-issue44-readback.json').write_text(json.dumps(readback,indent=2)+'\n')
comment='Scoped post-IR MCJIT CodeGenNone outcome: '+paragraph+'\nAll original gates and public source0109ab46 remain unchanged; issue44 and the ultimate goal remain open.'
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
