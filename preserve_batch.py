"""Recover the finite retired diagnostic batch before any WORK-only removal."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,os,shutil,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
ARCHIVE=ROOT/'attestations/native-vulkan/2053ef7/native-shader-clock-placement'
HEAD='2053ef79c90752d37b0075b8137cfb9b571b6820'
TAG='evidence/issue-44-native-shader-clock-placement-2053ef7'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def git(*args,env=None,input=None):return subprocess.check_output(['git',*args],env=env,input=input,text=True,timeout=30).strip()
def birth(pid):
    try:return Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[19]
    except (FileNotFoundError,ProcessLookupError):return None
def guard():
    assert WORK==ROOT/'out/native-shader-clock-placement' and not WORK.is_symlink()
    assert git('rev-parse','HEAD')==HEAD==git('rev-parse','@{u}') and not git('status','--porcelain')
    assert git('rev-list','--left-right','--count','main...origin/main')=='0\t0'
    assert git('worktree','list','--porcelain').count('worktree ')==1
    assert read(WORK/'independent-placement-review.json')['accepted']
    for folder in ('','endpoint','transport','endpoint-query-budget'):
        w=WORK/folder
        assert read(w/'retirement-run.json')['accepted'] and read(w/'native-retirement.json')['pass']
        r=read(w/'placement-query-bridge.json')
        assert r['linux_child_birth_absent'] and r['linux_child_group_absent']
        assert birth(r['linux_child_pid'])!=r['linux_child_start_ticks']
        if folder in ('','endpoint'):assert not r['accepted'] and read(w/'placement-query-emergency.json')['recorded_births_absent']
        else:assert r['accepted']
        for path,v in read(w/'plan.json')['input_seals'].items():assert seal(ROOT/path)==v
    for p in read(WORK/'plan.json')['providers']:assert seal(Path(p['linux_path']))=={k:p[k] for k in ('bytes','sha256')}
    assert read(WORK/'ci-2053-final.json')['headSha']==HEAD and read(WORK/'ci-2053-final.json')['conclusion']=='success'
mode=sys.argv[1];assert mode in ('archive','close','remove');guard()
if mode=='archive':
    assert not ARCHIVE.exists();ARCHIVE.mkdir(parents=True)
    env=dict(os.environ);env['GIT_INDEX_FILE']=str(WORK/'preservation.index')
    useful=sorted(p for p in WORK.rglob('*') if p.is_file() and (p.suffix in ('.py','.cs','.ps1') or p.name in ('placement-analysis.json','independent-placement-review.json','source-review.json')))
    blobs={}
    try:
        git('read-tree','--empty',env=env)
        for p in useful:
            assert not p.is_symlink();name=str(p.relative_to(WORK));blob=git('hash-object','-w',str(p))
            git('update-index','--add','--cacheinfo','100644',blob,name,env=env);blobs[name]=blob
        tree=git('write-tree',env=env)
        commit=git('commit-tree',tree,'-p',HEAD,input='Preserve reviewed native shader-clock placement query sources (#44)\n\nTransport and Endpoint static queries bind exact diagnostic modules and ordinary driver code; zero shader dispatches. Initial512MiB attempts invalidated; Endpoint query uses the pre-existing90s/4GiB compiler-tooling budget. No scientific/performance/release verdict.\n')
        git('tag',TAG,commit)
    finally:
        for p in WORK.glob('preservation.index*'):p.unlink()
    rows=[]
    for p in sorted(p for p in WORK.rglob('*') if p.is_file()):
        assert not p.is_symlink();dest=ARCHIVE/'files'/p.relative_to(WORK);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
        row={'origin':str(p.relative_to(ROOT)),'recovery':str(dest.relative_to(ROOT)),**seal(p)}
        if str(p.relative_to(WORK)) in blobs:row['git_path']=str(p.relative_to(WORK));row['git_blob']=blobs[row['git_path']]
        rows.append(row)
    manifest={'source_revision':HEAD,'source_tag':TAG,'source_commit':commit,'files':len(rows),'bytes':sum(r['bytes'] for r in rows),'records':rows,
              'scope':'Complete finiteWORK recovery; useful source and exact review/placement findings in Git. Raw binary/driver/ownership evidence remains in ignored local attestations; no remote raw-backup claim. Recorded/scoped ownership only; hidden/global/scheduled absence unclaimed.'}
    (ARCHIVE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:manifest[k] for k in ('source_revision','source_tag','source_commit','files','bytes')}))
elif mode=='close':
    assert read(WORK/'preservation-review.json')['accepted']
    assert read(WORK/'publication-readback.json')['exact']
    assert read(WORK/'final-retirement-run.json')['accepted'] and read(WORK/'final-native-retirement.json')['pass']
    manifest=read(ARCHIVE/'manifest.json');assert git('ls-remote','origin','refs/tags/'+TAG).split()[0]==manifest['source_commit']
    rows=[]
    for p in sorted(p for p in WORK.rglob('*') if p.is_file()):
        assert not p.is_symlink();dest=ARCHIVE/'files'/p.relative_to(WORK);dest.parent.mkdir(parents=True,exist_ok=True)
        if dest.exists():assert seal(dest)==seal(p)
        else:shutil.copyfile(p,dest)
        rows.append({'origin':str(p.relative_to(ROOT)),'recovery':str(dest.relative_to(ROOT)),**seal(p)})
    for row in manifest['records']:assert seal(ROOT/row['recovery'])=={k:row[k] for k in ('bytes','sha256')}
    pub=ARCHIVE/'publication';pub.mkdir()
    (pub/'closing-ledger.json').write_text(json.dumps({'manifest':seal(ARCHIVE/'manifest.json'),'source_tag':TAG,'source_commit':manifest['source_commit'],'records':rows},indent=2)+'\n')
    print(json.dumps({'files':len(rows),'bytes':sum(r['bytes'] for r in rows),'ledger':seal(pub/'closing-ledger.json')}))
else:
    ledger=read(ARCHIVE/'publication/closing-ledger.json')
    assert read(WORK/'preservation-review.json')['accepted'] and read(WORK/'publication-readback.json')['exact']
    assert read(WORK/'final-retirement-run.json')['accepted'] and read(WORK/'final-native-retirement.json')['pass']
    assert git('ls-remote','origin','refs/tags/'+TAG).split()[0]==ledger['source_commit']
    expected={r['origin']:{k:r[k] for k in ('bytes','sha256')} for r in ledger['records']}
    assert {str(p.relative_to(ROOT)):seal(p) for p in WORK.rglob('*') if p.is_file()}==expected
    for row in ledger['records']:assert seal(ROOT/row['recovery'])=={k:row[k] for k in ('bytes','sha256')}
    shutil.rmtree(WORK)
    state={'source_revision':git('rev-parse','HEAD'),'clean':not git('status','--porcelain'),
           'development_tracking':git('rev-list','--left-right','--count','HEAD...@{u}'),'main_tracking':git('rev-list','--left-right','--count','main...origin/main'),
           'worktrees':git('worktree','list'),'work_absent':not WORK.exists(),'out_entries':sorted(p.name for p in (ROOT/'out').iterdir()),
           'removed_files':len(expected),'removed_bytes':sum(v['bytes'] for v in expected.values()),'committed_and_pushed':True,
           'scope':'Finite diagnostic closure only; scientific rendering, CPU integration and release qualification remain open.'}
    assert state['clean'] and state['development_tracking']==state['main_tracking']=='0\t0'
    (ARCHIVE/'publication/final-state.json').write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state))
