"""Preserve this completed trial's tooling without changing development HEAD."""
from pathlib import Path
import hashlib,json,os,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
source=json.loads((WORK/'source.json').read_text());head=source['source_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'])
tag='evidence/issue-90-restoration-tools-'+head[:7];index=WORK/'source-git.index'
assert not index.exists()
assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+tag]).returncode==1
paths=[p for p in sorted(WORK.rglob('*')) if p.is_file() and p.suffix in ('.py','.ps1') and '__pycache__' not in p.parts and 'temp' not in p.relative_to(WORK).parts]
assert paths
records=[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths]
env=dict(os.environ,GIT_INDEX_FILE=str(index))
def git(*args,**kwargs):return subprocess.check_output(['git',*args],env=env,**kwargs)
try:
 git('read-tree',head);git('add','-f','--',*[r['path'] for r in records])
 tree=git('write-tree',text=True).strip()
 commit=git('commit-tree',tree,'-p',head,input='Preserve rejected paired Transport restoration and terminal tools (#90)\n',text=True).strip()
 for r in records:
  data=git('show',commit+':'+r['path']);assert len(data)==r['bytes'] and hashlib.sha256(data).hexdigest()==r['sha256']
 git('tag',tag,commit);git('push','origin','refs/tags/'+tag)
 assert subprocess.check_output(['git','ls-remote','origin','refs/tags/'+tag],text=True).split()[0]==commit
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
 assert not subprocess.check_output(['git','status','--porcelain'])
 receipt={'tag':tag,'commit':commit,'base':head,'sources':records,'source_count':len(records),'pushed_verified':True,'scope':'Source-only evidence preservation; development HEAD unchanged. Actual gate disposition recorded separately.'}
 q=WORK/'source-git-preservation.json';assert not q.exists();q.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='sources'}))
finally:index.unlink(missing_ok=True)
