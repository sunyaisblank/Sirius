from pathlib import Path
import os,subprocess,json,hashlib
r=Path.cwd();w=Path(__file__).resolve().parent
head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();assert head=='d92310c18ccb6141cf5c967f61be706f4603edd9' and not subprocess.check_output(['git','status','--porcelain'])
excluded={'executed-inputs','whole-sealed-inputs','candidate-source','restored-source'}
paths=[p for p in sorted(w.rglob('*.py')) if not excluded.intersection(p.relative_to(w).parts)]
assert paths
index=w/'source-preservation-index';assert not index.exists();env=os.environ.copy();env['GIT_INDEX_FILE']=str(index)
def run(*args):return subprocess.check_output(['git',*args],env=env,text=True).strip()
try:
 run('read-tree',head);records={}
 for p in paths:
  blob=run('hash-object','-w',str(p));relative=str(p.relative_to(r));run('update-index','--add','--cacheinfo','100644,'+blob+','+relative);records[relative]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'git_blob':blob}
 tree=run('write-tree');commit=subprocess.check_output(['git','commit-tree',tree,'-p',head],input='Preserve exact Dense trial and restoration controls (#88)\n',env=env,text=True).strip();tag='evidence/issue-88-tools-d92310c';run('update-ref','refs/tags/'+tag,commit)
 for name,seal in records.items():assert hashlib.sha256(subprocess.check_output(['git','show',commit+':'+name])).hexdigest()==seal['sha256']
 subprocess.run(['git','push','origin','refs/tags/'+tag],check=True)
 remote=subprocess.check_output(['git','ls-remote','--tags','origin','refs/tags/'+tag],text=True).strip().split()[0];assert remote==commit
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head and not subprocess.check_output(['git','status','--porcelain'])
 report={'source_commit':commit,'source_tree':tree,'parent_production_revision':head,'tag':tag,'remote_exact':True,'paths':records,'unique_bodies':len({q['sha256'] for q in records.values()}),'production_head_unchanged':True,'scope':'Exact useful task-generated source tools preserved in pushed historical source-only Git snapshot. Production/source files already preserved in their commits. No raw/binary evidence remote-backup claim.'};(w/'source-preservation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['source_commit','tag','unique_bodies','remote_exact']}))
finally:
 index.unlink(missing_ok=True);Path(str(index)+'.lock').unlink(missing_ok=True)
