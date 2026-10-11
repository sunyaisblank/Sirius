"""Preserve useful diagnostic source and raw receipts; remove only the finished workspace."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, os, shutil, subprocess
ROOT = Path.cwd(); WORK = Path(__file__).resolve().parent
HEAD = '451e2cbf833b87f51b1819421d68163a4b5fd1e5'
BASELINE = '2053ef79c90752d37b0075b8137cfb9b571b6820'
ARCHIVE = ROOT/'attestations/native-vulkan/451e2cb/native-division-phase-observer'
NATIVE = ROOT/'bin/windows-msvc/native-451e2cb-division-phase-observer'
TAG = 'evidence/issue-44-native-division-phase-observer-451e2cb'
PREFIX = 'diagnostics/native-division-phase-observer/'
CASES = ('device-timestamps','mixed-transport','factored-endpoints','placement')

def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):
    with p.open('rb') as stream: digest = hashlib.file_digest(stream,'sha256').hexdigest()
    return {'bytes':p.stat().st_size,'sha256':digest}
def git(*args,env=None,input=None):
    return subprocess.check_output(['git',*args],env=env,input=input,text=True,timeout=30).strip()
def files(base):
    entries = list(base.rglob('*'))
    assert all(not p.is_symlink() and (p.is_file() or p.is_dir()) for p in entries)
    return sorted(p for p in entries if p.is_file())
def copy_all():
    records = []
    for p in files(WORK):
        target = ARCHIVE/'files'/p.relative_to(WORK);target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists(): assert seal(target) == seal(p)
        else: shutil.copyfile(p,target)
        assert seal(p) == seal(target)
        records.append({'origin':str(p.relative_to(ROOT)),'recovery':str(target.relative_to(ROOT)),**seal(p)})
    return records

def common_guard():
    assert ROOT.resolve() == Path('/home/astra/.project/Sirius')
    assert WORK == ROOT/'out/native-division-phase-observer' and not WORK.is_symlink()
    assert not git('status','--porcelain')
    assert git('rev-list','--left-right','--count','HEAD...@{u}') == '0\t0'
    assert git('rev-list','--left-right','--count','main...origin/main') == '0\t0'
    assert git('worktree','list','--porcelain').count('worktree ') == 1
    assert git('rev-parse','development/retained-dense-step-control-2026-10-06') == BASELINE
    assert git('rev-list','--left-right','--count','development/retained-dense-step-control-2026-10-06...origin/development/retained-dense-step-control-2026-10-06') == '0\t0'
    assert read(WORK/'actual-results-review.json')['accepted']
    assert read(WORK/'results.json')['accepted']
    review = read(WORK/'preservation-source-review.json')
    assert review['accepted'] and review['helper_seal'] == seal(Path(__file__))
    for case in CASES:
        directory = WORK/case
        bridge = read(directory/('placement-query-bridge.json' if case == 'placement' else 'native-control-bridge.json'))
        assert bridge['accepted'] and read(directory/'retirement-run.json')['accepted'] and read(directory/'native-retirement.json')['pass']
    ci = read(WORK/'ci-completed.json')
    assert ci['headSha'] == HEAD and ci['conclusion'] == 'success'
    assert read(WORK/'actual-placement-review.json')['accepted']
    assert read(WORK/'actual-source-build-review.json')['accepted']

assert len(sys.argv) == 2 and sys.argv[1] in ('archive','remove')
common_guard()
if sys.argv[1] == 'archive':
    assert git('rev-parse','HEAD') == HEAD
    assert not ARCHIVE.exists() and not git('tag','--list',TAG)
    plan = read(WORK/'execution-plan.json');assert plan['source_revision'] == HEAD
    for path,expected in plan['input_seals'].items(): assert seal(ROOT/path) == expected
    for path,expected in read(WORK/'actual-results-review.json')['reviewed_seals'].items(): assert seal(Path(path)) == expected
    env = dict(os.environ,GIT_INDEX_FILE=str(WORK/'preservation.index'));blobs = {}
    try:
        git('read-tree',HEAD,env=env)
        for p in files(WORK):
            if p.suffix not in ('.py','.ps1','.cs','.slang'): continue
            if {'owner-temp','native-temp'}.intersection(p.relative_to(WORK).parts): continue
            relative = str(p.relative_to(WORK));blob = git('hash-object','-w',str(p));name = PREFIX+relative
            git('update-index','--add','--cacheinfo','100644',blob,name,env=env)
            blobs[relative] = {'git_path':name,'git_blob':blob}
        tree = git('write-tree',env=env)
        message = (WORK/'preservation-commit-message.txt').read_text()
        assert message.startswith('Preserve native division observations and source (#44)\n')
        commit = git('commit-tree',tree,'-p',HEAD,input=message)
        git('tag',TAG,commit)
    finally:
        for p in WORK.glob('preservation.index*'): p.unlink()
    ARCHIVE.mkdir(parents=True);records = copy_all()
    for row in records: row.update(blobs.get(str(Path(row['origin']).relative_to(WORK.relative_to(ROOT))),{}))
    manifest = {'source_revision':HEAD,'restoration_revision':BASELINE,'source_tag':TAG,'source_commit':commit,'records':records,
        'scope':'Complete diagnostic source plus useful preparers/controllers/readers in Git; concise findings in commit and issue44. Full raw receipts/artifacts retained locally in ignored attestations. No remote raw backup or full qualification claim.'}
    (ARCHIVE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'source_tag':TAG,'source_commit':commit,'files':len(records),'bytes':sum(r['bytes'] for r in records)}))
else:
    assert git('rev-parse','HEAD') == BASELINE and git('branch','--show-current') == 'development/retained-dense-step-control-2026-10-06'
    review = read(WORK/'cleanup-ready-review.json')
    assert review['accepted'] and review['helper_seal'] == seal(Path(__file__))
    assert read(WORK/'publication-readback.json')['exact']
    assert read(WORK/'restored-source-build.json')['accepted'] and read(WORK/'restored-products.json')['accepted']
    bridge = read(WORK/'final-retirement-run.json');census = read(WORK/'final-native-retirement.json')
    assert bridge['accepted'] and bridge['inputs_unchanged'] and not bridge['errors']
    assert census['pass'] and not census['matches'] and not census['errors']
    assert len(census['targets']) == 8 and all(t['recorded_birth_absent'] for t in census['targets'])
    for path,expected in bridge['input_seals'].items(): assert seal(Path(path)) == expected
    assert census['native_stage'] == subprocess.check_output(['wslpath','-w',str(NATIVE)],text=True,timeout=10).strip()
    assert census['stage'] == subprocess.check_output(['wslpath','-w',str(WORK)],text=True,timeout=10).strip()
    targets = []
    for case in CASES:
        owner = read(WORK/case/'run/owner.json')
        targets.extend([(case,'owner',int(owner['owner_pid']),owner['owner_start_time_utc']),
                        (case,'child',int(owner['test_pid']),owner['test_start_time_utc'])])
    assert len(set(targets)) == 8
    assert [(t['case'],t['kind'],int(t['pid']),t['birth']) for t in census['targets']] == targets
    manifest = read(ARCHIVE/'manifest.json')
    assert git('ls-remote','origin','refs/tags/'+TAG).split()[0] == manifest['source_commit']
    for row in manifest['records']: assert seal(ROOT/row['origin']) == seal(ROOT/row['recovery']) == {k:row[k] for k in ('bytes','sha256')}
    plan = read(WORK/'execution-plan.json');prefix = str(NATIVE.relative_to(ROOT))+'/'
    expected = {p:v for p,v in plan['input_seals'].items() if p.startswith(prefix)}
    assert {str(p.relative_to(ROOT)):seal(p) for p in files(NATIVE)} == expected
    retained = {(row['bytes'],row['sha256']) for row in manifest['records']}
    assert all((v['bytes'],v['sha256']) in retained for v in expected.values()), 'Stage file has no archived counterpart'
    births = []
    def absent(pid,ticks):
        try: actual = Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[19]
        except (FileNotFoundError,ProcessLookupError): actual = None
        assert ticks and actual != ticks
        births.append({'pid':pid,'start_ticks':ticks,'recorded_birth_absent':True})
    for case in CASES:
        directory = WORK/case
        b = read(directory/('placement-query-bridge.json' if case == 'placement' else 'native-control-bridge.json'))
        rr = read(directory/'retirement-run.json')
        absent(b['linux_owner_pid'],b['linux_owner_start_ticks']);absent(b['linux_child_pid'],b['linux_child_start_ticks']);absent(rr['pid'],rr['start_ticks'])
    matches = [];hidden = [];visible = 0
    for p in Path('/proc').iterdir():
        if not p.name.isdigit() or int(p.name) in (os.getpid(),os.getppid()): continue
        try: command = (p/'cmdline').read_bytes().replace(b'\0',b' ');visible += 1
        except (FileNotFoundError,ProcessLookupError,PermissionError) as error:
            hidden.append({'pid':int(p.name),'reason':type(error).__name__});continue
        if any(v.encode() in command for v in (str(WORK),str(WORK.relative_to(ROOT)),str(NATIVE),str(NATIVE.relative_to(ROOT)))):
            matches.append({'pid':int(p.name),'command':command.decode(errors='replace')})
    assert not matches,matches
    records = copy_all();publication = ARCHIVE/'publication';publication.mkdir()
    ledger = {'source_tag':TAG,'source_commit':manifest['source_commit'],'manifest_seal':seal(ARCHIVE/'manifest.json'),'records':records,
        'native_stage_files':len(expected),'native_stage_bytes':sum(v['bytes'] for v in expected.values()),'linux_births':births,
        'linux_visible_processes':visible,'linux_hidden_records':hidden,
        'scope':'Exact saved births and fresh accessible references only. No global/hidden/scheduled absence claim.'}
    (publication/'closing-ledger.json').write_text(json.dumps(ledger,indent=2)+'\n')
    shutil.rmtree(NATIVE);shutil.rmtree(WORK)
    state = {'source_revision':BASELINE,'clean':not git('status','--porcelain'),'tracking':git('rev-list','--left-right','--count','HEAD...@{u}'),
        'main_tracking':git('rev-list','--left-right','--count','main...origin/main'),'worktrees':git('worktree','list'),
        'finished_work_absent':not WORK.exists(),'finished_native_stage_absent':not NATIVE.exists(),
        'out_entries':sorted(p.name for p in (ROOT/'out').iterdir()),'removed_work_files':len(records),'removed_work_bytes':sum(v['bytes'] for v in records),
        'scope':'Finite diagnostic source/raw evidence preserved; accepted source/products restored. Long-term physical rendering and science/release qualification remain open.'}
    assert state['clean'] and state['tracking'] == state['main_tracking'] == '0\t0'
    (publication/'removal-state.json').write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state))
