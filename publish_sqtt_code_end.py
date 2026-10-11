"""Publish the reviewed finite findings to the existing GitHub remediation issue."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,subprocess
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def gh(*a):return subprocess.check_output(['gh',*a],text=True).strip()
source=read(W/'source-preservation.json');assert source['pushed'] and source['remote_commit_verified']
assert subprocess.check_output(['git','ls-remote','origin','refs/tags/'+source['tag']],text=True).split()[0]==source['commit']
protection=read(W/'protection.json');review=read(W/'preservation-manifest-review.json');assert review['accepted'] and review['manifest']==protection['manifest']
before=read(W/'issue-44-before.json');current=json.loads(gh('issue','view','44','--json','number,state,body,url'));assert current['number']==44 and current['state']=='OPEN' and current['body']==before['body']
start=current['body'].index('**Next bounded investigation.**');end=current['body'].index('**Finite closure criteria.**',start)
replacement=(W/'next-investigation.md').read_text()
body=current['body'][:start]+replacement+current['body'][end:];p=W/'issue-44-body.md';assert not p.exists();p.write_text(body)
comment=W/'issue-44-update.md';assert not (W/'publication.json').exists()
inputs={str(q.relative_to(ROOT)):seal(q) for q in (Path(__file__),comment,W/'next-investigation.md',W/'issue-44-before.json',W/'source-preservation.json',W/'protection.json',W/'preservation-manifest-review.json')}
url=gh('issue','comment','44','--body-file',str(comment));comment_id=url.rsplit('-',1)[-1];assert comment_id.isdigit()
actual=json.loads(gh('api','repos/sunyaisblank/Sirius/issues/comments/'+comment_id));assert actual['body']==comment.read_text() and actual['html_url']==url
gh('issue','edit','44','--body-file',str(p));after=json.loads(gh('issue','view','44','--json','number,state,body,url'));assert after['body']==body and after['state']=='OPEN'
(W/'issue-44-after.json').write_text(json.dumps(after,indent=2)+'\n')
(W/'publication.json').write_text(json.dumps({'input_seals':inputs,'source_snapshot':source,'comment_url':url,'comment_readback_exact':True,'issue_body_readback_exact':True,'issue_open':True,'comment_seal':seal(comment),'body_seal':seal(p),'scope':'Existing issue44 reviewed findings and next bounded prerequisite only; no issue closure, merge or full-qualification claim.'},indent=2)+'\n');print(json.dumps({'comment_url':url,'comment_readback_exact':True,'body_readback_exact':True,'issue_open':True}))
