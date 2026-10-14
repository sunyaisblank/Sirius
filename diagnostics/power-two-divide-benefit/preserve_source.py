"""Preserve the frozen original benefit apparatus without changing source HEAD."""
from pathlib import Path
import hashlib,json,os,subprocess

R=Path.cwd().resolve();Q=Path(__file__).resolve().parent;W=Q/'benefit'
revision='4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
tag='evidence/issue-99-power-two-divide-benefit-4bcd987'
query=lambda args,**kwargs:subprocess.check_output(args,**kwargs)
assert query(['git','rev-parse','HEAD'],text=True).strip()==revision
assert not query(['git','status','--porcelain'],text=True)
assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+tag]).returncode==1
index=Q/'source.index';assert not index.exists()
environment=dict(os.environ,GIT_INDEX_FILE=str(index))
sources=sorted(Q.glob('*.py'))+sorted(W.glob('*.py'))
sources += [Q/'audited-observer.cpp']
sources += [p for p in sorted(Q.glob('*.json')) if not p.name.startswith('issue99-')]
sources += [p for p in sorted(W.glob('*.json')) if p.name != 'software-sequence.json']
assert (Q/'trial-acceptance.json').is_file() and (Q/'disposition.json').is_file()
records=[]
try:
    subprocess.run(['git','read-tree',revision],env=environment,check=True)
    for path in sources:
        relative='diagnostics/power-two-divide-benefit/'+str(path.relative_to(Q))
        blob=query(['git','hash-object','-w',str(path)],text=True).strip()
        subprocess.run(['git','update-index','--add','--cacheinfo','100644,'+blob+','+relative],env=environment,check=True)
        records.append({'path':relative,'origin':str(path.relative_to(R)),'blob':blob,'bytes':path.stat().st_size,
                        'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    tree=query(['git','write-tree'],env=environment,text=True).strip()
    commit=query(['git','commit-tree',tree,'-p',revision],input='Preserve issue #99 unchanged original matched benefit trial apparatus\n',text=True).strip()
    subprocess.run(['git','update-ref','refs/tags/'+tag,commit],check=True)
finally:index.unlink(missing_ok=True)
subprocess.run(['git','push','origin','refs/tags/'+tag],check=True)
remote=query(['git','ls-remote','origin','refs/tags/'+tag],text=True).strip()
assert remote==commit+'\trefs/tags/'+tag
assert query(['git','rev-parse','HEAD'],text=True).strip()==revision
assert not query(['git','status','--porcelain'],text=True)
(Q/'preserved-source.json').write_text(json.dumps({'tag':tag,'commit':commit,'parent':revision,'source_records':records,
    'remote_readback':remote,'working_tree_and_head_unchanged':True},indent=2)+'\n')
print(json.dumps({'tag':tag,'commit':commit,'source_files':len(records),'pushed':True}))
