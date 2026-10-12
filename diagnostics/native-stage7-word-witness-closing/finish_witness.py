"""Preserve closing source, recover consumed versions, and remove the finite witness.

Deletion requires independent clearance, pushed immutable source, exact local
recovery, restored products, terminal owners and fresh accessible process checks.
"""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, os, shutil, subprocess

ROOT = Path('/home/astra/.project/Sirius')
WORK = ROOT/'out/native-stage-seven-value-witness'
ARCHIVE = ROOT/'attestations/native-vulkan/5f27ddc/stage-seven-word-witness'
HEAD = '9fbd0e8b0d549a13476d066f2877dcd2cf9a94d9'
TRIAL = '5f27ddcfd8cc776a18dddf530e49cf905a2646b2'
NATIVE = ROOT/'bin/windows-msvc/native-5f27ddc-stage-seven-word-witness'
BRANCH = 'diagnostic/stage-seven-word-witness-2026-10-12'
TAG = 'evidence/issue-44-native-stage7-word-witness-closing-9fbd0e8-v2'

def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):
    with p.open('rb') as stream: digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': p.stat().st_size, 'sha256': digest}
def git(*args, **kw): return subprocess.check_output(['git', *args], cwd=ROOT, timeout=60, **kw)
def save(p, value):
    assert not p.exists(), p
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(value, indent=2)+'\n')
def files(base):
    assert base.is_dir() and not base.is_symlink() and base.resolve() == base
    paths = sorted(base.rglob('*'))
    assert all(not p.is_symlink() and (p.is_file() or p.is_dir()) for p in paths)
    return [p for p in paths if p.is_file()]
def state():
    assert Path.cwd().resolve() == ROOT and Path(__file__).resolve() == WORK/'finish_witness.py'
    assert git('rev-parse', 'HEAD').decode().strip() == HEAD
    assert git('branch', '--show-current').decode().strip() == 'development/retained-dense-step-control-2026-10-06'
    assert not git('status', '--porcelain')
    for refs in ('HEAD...@{u}', 'main...origin/main'):
        assert git('rev-list', '--left-right', '--count', refs).split() == [b'0', b'0']
    assert git('worktree', 'list', '--porcelain').decode().count('worktree ') == 1
def pushed(record):
    assert record['accepted'] and record['pushed']
    assert git('ls-remote', 'origin', 'refs/tags/'+record['source_tag']).split()[0].decode() == record['source_commit']
    for name, version in record['source_versions'].items():
        assert git('show', record['source_commit']+':'+version['git_path']) == (WORK/name).read_bytes()
def restored():
    r = read(WORK/'restoration-products.json')
    assert r['source_revision'] == HEAD
    for name, expected in r['products'].items():
        assert expected['matches_saved_baseline'] and seal(ROOT/name) == {k:expected[k] for k in ('bytes','sha256')}
    assert seal(ROOT/'bin/linux-gcc/src/sirius/app/sirius_app_tests') == r['app_tests']
    receipt = read(ROOT/'bin/linux-gcc/generated/sirius/alignment_receipt.json')
    assert receipt == r['generated_receipt'] and len(receipt['admitted_domains']) == 2
    assert len(receipt['pending_domains']) == 6 and not receipt['ultimate_ideal']['satisfied']
def versions():
    """Join all recorded owned inputs; archives include historical aliases."""
    wanted = {}
    ledgers = []
    for p in sorted(WORK.rglob('inputs-before.json')):
        before = read(p); assert before == read(p.with_name('inputs-after.json'))
        ledgers.append((str(p.relative_to(WORK)), before))
    for p in sorted(WORK.rglob('*.json')):
        if p.name.startswith('inputs-') or p.name == 'input-recovery.json': continue
        r = read(p)
        if isinstance(r, dict) and isinstance(r.get('input_seals'), dict):
            ledgers.append((str(p.relative_to(WORK)), r['input_seals']))
    for ledger, records in ledgers:
        for name, value in records.items():
            if not isinstance(value, dict) or not {'bytes','sha256'} <= value.keys(): continue
            key = (name, value['bytes'], value['sha256'])
            wanted.setdefault(key, []).append(ledger)
    local = {}
    for p in files(ARCHIVE):
        s = seal(p); local.setdefault((s['bytes'],s['sha256']), str(p.relative_to(ROOT)))
    tracked = set(git('ls-tree','-r','--name-only',TRIAL).decode().splitlines())
    source_paths={}
    for filename in ('source-preservation-pushed.json','closing-source-preservation.json','closing-source-preservation-v2.json'):
        if not (WORK/filename).exists(): continue
        source=read(WORK/filename)
        for name,r in source['source_versions'].items():
            origin=str((WORK/name).relative_to(ROOT))
            source_paths.setdefault(origin,[]).append((source['source_commit'],r['git_path']))
    result = []
    for (name, size, digest), ledgers in sorted(wanted.items()):
        recovery = None
        rel = str(Path(name).relative_to(ROOT)) if Path(name).is_absolute() and Path(name).is_relative_to(ROOT) else name
        if rel in tracked:
            for revision in (TRIAL, HEAD):
                data = git('show',revision+':'+rel)
                if len(data) == size and hashlib.sha256(data).hexdigest() == digest:
                    recovery = 'git:'+revision+':'+rel; break
        if recovery is None:
            for revision,path in source_paths.get(rel,[]):
                data=git('show',revision+':'+path)
                if len(data)==size and hashlib.sha256(data).hexdigest()==digest:
                    recovery='git:'+revision+':'+path;break
        if recovery is None: recovery = local.get((size,digest))
        if recovery is None:
            p = Path(name) if Path(name).is_absolute() else ROOT/name
            assert p.exists() and not p.is_symlink() and not p.is_relative_to(WORK) and not p.is_relative_to(NATIVE), name
            assert seal(p) == {'bytes':size,'sha256':digest}, name
            recovery = 'retained:'+str(p)
        result.append({'origin':name,'bytes':size,'sha256':digest,'recovery':recovery,'ledgers':ledgers})
    return result
def birth(pid):
    try:
        fields = Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()
        return fields[19], int(fields[1])
    except (FileNotFoundError,ProcessLookupError): return None,None
def terminal_births():
    boot = Path('/proc/sys/kernel/random/boot_id').read_text().strip(); saved = set()
    for p in WORK.rglob('owner.json'):
        r = read(p)
        if 'boot_id' not in r: continue
        assert r['boot_id'] == boot and r['passed'] and r['status'] == 'terminal'
        assert r['child_exit'] == 0 and r['source_unchanged'] and r['whole_inputs_unchanged']
        assert r['stop_reason'] is None and not r['cleanup_errors'] and not r['remaining_owned_processes']
        saved.add((int(r['controller']['pid']),str(r['controller']['start_ticks'])))
        saved.update((int(t['pid']),str(t['saved_start_ticks'])) for t in r['observed_births'])
    for p in WORK.rglob('*.json'):
        r = read(p)
        if not isinstance(r,dict): continue
        if 'linux_owner_pid' in r:
            assert r['accepted']; saved.add((int(r['linux_owner_pid']),str(r['linux_owner_start_ticks'])))
        if r.get('accepted') and 'pid' in r and 'start_ticks' in r:
            assert r.get('birth_absent') and r.get('group_absent')
            saved.add((int(r['pid']),str(r['start_ticks'])))
        if r.get('accepted') and 'linux_child_pid' in r:
            saved.add((int(r['linux_child_pid']),str(r['linux_child_start_ticks'])))
    result = []
    for pid,ticks in sorted(saved):
        current,_ = birth(pid); assert current != ticks,(pid,ticks)
        result.append({'pid':pid,'saved_start_ticks':ticks,'current_start_ticks':current,'recorded_birth_absent':True})
    return result
def consumers():
    ancestors=set(); pid=os.getpid()
    while pid and pid not in ancestors:
        ancestors.add(pid); _,pid=birth(pid)
    visible=[]; inaccessible=[]; matches=[]
    terms=[str(p).encode() for base in (WORK,NATIVE) for p in (base,base.relative_to(ROOT))]
    for proc in sorted(Path('/proc').iterdir()):
        if not proc.name.isdecimal() or int(proc.name) in ancestors: continue
        try:
            command=(proc/'cmdline').read_bytes(); cwd=(proc/'cwd').resolve(strict=True)
        except (PermissionError,FileNotFoundError,ProcessLookupError) as e:
            inaccessible.append({'pid':int(proc.name),'reason':type(e).__name__}); continue
        visible.append(int(proc.name))
        if any(t in command for t in terms) or any(cwd==b or b in cwd.parents for b in (WORK,NATIVE)):
            matches.append(int(proc.name))
    assert not matches,matches
    return {'accessible_pids':visible,'inaccessible_or_disappeared':inaccessible,'excluded_ancestors':sorted(ancestors),'matches':matches}

assert len(sys.argv)==2 and sys.argv[1] in ('source','recovery','remove')
state(); restored()
assert read(WORK/'actual-results-review.json')['accepted']
assert read(WORK/'final-retirement-run.json')['accepted'] and read(WORK/'final-native-retirement.json')['pass']
initial=read(WORK/'source-preservation-pushed.json'); pushed(initial)
if sys.argv[1]=='source':
    assert not git('tag','--list',TAG).strip()
    index=WORK/'closing-source.index'; assert not index.exists()
    env=dict(os.environ,GIT_INDEX_FILE=str(index)); records={}
    try:
        parent=read(WORK/'closing-source-preservation.json')['source_commit']
        git('read-tree',parent,env=env)
        for name in ('prepare_final_retirement.py','run_final_retirement.py','verify_final_retirement.ps1','run_cleanup_retirement.py','finish_witness.py','versions/finish_witness.v1.py','versions/finish_witness.v2.py'):
            p=WORK/name; target='diagnostics/native-stage7-word-witness-closing/'+name
            blob=git('hash-object','-w','--stdin',input=p.read_bytes()).decode().strip()
            git('update-index','--add','--cacheinfo','100644,'+blob+','+target,env=env)
            records[name]={'git_path':target,**seal(p)}
        tree=git('write-tree',env=env).decode().strip()
        message='Preserve completed native word witness closing source (#44)\n\nAll three original controls pass;162 curved high pairs match5022 complete five-word results. Production9fbd and six baseline products are restored byte-exact; qualification remains2/8 build-only. Preserve bounded retirement and exact input recovery/deletion guards. No optimization, rendering gain, universal word theorem or full qualification is admitted.\n'
        commit=git('commit-tree',tree,'-p',parent,input=message.encode()).decode().strip()
        git('update-ref','refs/tags/'+TAG,commit); git('push','origin','refs/tags/'+TAG)
        r={'accepted':True,'pushed':True,'source_revision':HEAD,'source_tag':TAG,'source_commit':commit,'source_versions':records}
        pushed(r);save(WORK/'closing-source-preservation-v2.json',r)
    finally: index.unlink(missing_ok=True)
    print(json.dumps({'source_tag':TAG,'source_commit':commit,'useful_sources':len(records)}))
elif sys.argv[1]=='recovery':
    pushed(read(WORK/'closing-source-preservation-v2.json'))
    recovered=versions()
    save(WORK/'input-recovery.json',{'accepted':True,'records':recovered,'origin_versions':len(recovered),'missing':0,'scope':'All recorded owned input seals, including historical aliases, diagnostic products and restored current inputs.'})
    print(json.dumps({'recovered_origin_versions':len(recovered)}))
else:
    clearance=read(WORK/'cleanup-clearance.json'); assert clearance['accepted']
    for name,expected in clearance['approved_seals'].items(): assert seal(WORK/name)==expected
    assert clearance['approved_seals']['finish_witness.py']==seal(Path(__file__))
    closing=read(WORK/'closing-source-preservation-v2.json'); pushed(closing)
    assert read(WORK/'publication-accepted.json')['accepted']
    for name,expected in read(WORK/'publication-accepted.json')['readback_seals'].items(): assert seal(WORK/name)==expected
    assert read(WORK/'input-recovery.json')['records']==versions()
    bridge=read(WORK/'final-retirement-run.json'); census=read(WORK/'final-native-retirement.json')
    assert bridge['inputs_unchanged'] and not bridge['errors'] and not bridge['signals']
    assert not census['matches'] and not census['errors'] and len(census['targets'])==6
    assert all(t['recorded_birth_absent'] for t in census['targets'])
    for name,expected in bridge['input_seals'].items(): assert seal(Path(name))==expected
    native_expected={n:v for n,v in read(WORK/'execution-plan.json')['input_seals'].items() if n.startswith(str(NATIVE.relative_to(ROOT))+'/')}
    assert len(native_expected)==83 and {str(p.relative_to(ROOT)):seal(p) for p in files(NATIVE)}==native_expected
    for refs in (BRANCH+'...origin/'+BRANCH,): assert git('rev-list','--left-right','--count',refs).split()==[b'0',b'0']
    assert git('rev-parse',BRANCH).decode().strip()==TRIAL
    records=[]
    directory_sets={str(base.relative_to(ROOT)):sorted(str(p.relative_to(base)) for p in base.rglob('*') if p.is_dir()) for base in (WORK,NATIVE)}
    for base,part in ((WORK,'files'),(NATIVE,'native-stage')):
        for p in files(base):
            q=ARCHIVE/part/p.relative_to(base); expected=seal(p)
            if p==Path(__file__).resolve(): q=ARCHIVE/'closing-v2/finish_witness.py'
            if q.exists(): assert seal(q)==expected,(p,q)
            else: q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
            assert seal(p)==seal(q)==expected
            records.append({'origin':str(p.relative_to(ROOT)),'recovery':str(q.relative_to(ROOT)),**expected})
    # Fresh observations follow expensive recovery/copying. The wrapper retains
    # the original 40s ownership/retirement deadline and sealed PS source.
    assert not (WORK/'cleanup-retirement-run.json').exists()
    with (WORK/'cleanup-retirement-controller.stdout').open('xb') as out,(WORK/'cleanup-retirement-controller.stderr').open('xb') as err:
        subprocess.check_call([sys.executable,'-B',str(WORK/'run_cleanup_retirement.py')],stdout=out,stderr=err,cwd=ROOT)
    fresh=read(WORK/'cleanup-retirement-run.json');fresh_census=read(WORK/'cleanup-native-retirement.json')
    assert fresh['accepted'] and fresh['inputs_unchanged'] and not fresh['errors'] and not fresh['signals']
    assert fresh_census['pass'] and not fresh_census['matches'] and not fresh_census['errors']
    assert fresh_census['targets']==census['targets']
    archived={r['origin']:r for r in records}
    for name,expected in fresh['input_seals'].items():
        assert seal(Path(name))==expected
        recovery=archived[str(Path(name).relative_to(ROOT))]
        assert seal(ROOT/recovery['recovery'])==expected
    known={r['origin'] for r in records}
    for p in files(WORK):
        origin=str(p.relative_to(ROOT))
        if origin in known: continue
        q=ARCHIVE/'files'/p.relative_to(WORK);expected=seal(p)
        assert not q.exists();q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
        assert seal(q)==expected
        records.append({'origin':origin,'recovery':str(q.relative_to(ROOT)),**expected})
    # Refuse changed bytes, new entries or root/symlink replacement between
    # copying and removal. Fresh observer outputs have already joined records.
    for base in (WORK,NATIVE):
        expected={r['origin']:{k:r[k] for k in ('bytes','sha256')} for r in records if (ROOT/r['origin']).is_relative_to(base)}
        assert {str(p.relative_to(ROOT)):seal(p) for p in files(base)}==expected
        assert sorted(str(p.relative_to(base)) for p in base.rglob('*') if p.is_dir())==directory_sets[str(base.relative_to(ROOT))]
    state();restored()
    births=terminal_births();census_linux=consumers()
    save(ARCHIVE/'closing-ledger.json',{'accepted':True,'records':records,'directories':directory_sets,'linux_births':births,'linux_consumers':census_linux,'source_tag':closing['source_tag'],'source_commit':closing['source_commit'],'scope':'Complete local byte recovery and pushed useful source; finite saved births/accessibility only. No hidden/global/scheduled-consumer absence claim.'})
    shutil.rmtree(NATIVE); shutil.rmtree(WORK)
    git('push','origin','--delete',BRANCH)
    git('branch','-D',BRANCH)
    assert not git('branch','--list',BRANCH).strip() and not git('ls-remote','origin','refs/heads/'+BRANCH).strip()
    assert not git('status','--porcelain')
    for refs in ('HEAD...@{u}','main...origin/main'): assert git('rev-list','--left-right','--count',refs).split()==[b'0',b'0']
    save(ARCHIVE/'cleanup-result.json',{'accepted':True,'source_revision':HEAD,'removed':[str(WORK.relative_to(ROOT)),str(NATIVE.relative_to(ROOT)),BRANCH],'preserved_files':len(records),'saved_linux_births':len(births),'out_entries':sorted(p.name for p in (ROOT/'out').iterdir()),'worktrees':git('worktree','list').decode(),'clean_and_pushed':True,'full_goal_complete':False,'qualification':'2/8 build-only; scientific/render/runtime/release obligations remain open.'})
    print(json.dumps({'removed_finished_witness':True,'preserved_files':len(records),'saved_linux_births':len(births)}))
