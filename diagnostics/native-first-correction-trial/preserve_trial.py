"""Preserve useful source and all raw trial versions before removing finished workspace."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, os, shutil, subprocess
ROOT = Path.cwd(); WORK = Path(__file__).resolve().parent; RUNTIME = WORK/'runtime'
TRIAL = 'ea6bbb03f8a7fe6fa9e3939342b058e25514c129'
BASELINE = '2053ef79c90752d37b0075b8137cfb9b571b6820'
BRANCH = 'development/retained-dense-step-control-2026-10-06'
ARCHIVE = ROOT/'attestations/native-vulkan/ea6bbb0/native-division-first-correction-trial'
NATIVE = ROOT/'bin/windows-msvc/native-ea6bbb0-division-first-correction-trial'
TAG = 'evidence/issue-44-native-first-correction-trial-ea6bbb0'
PREFIX = 'diagnostics/native-first-correction-trial/'
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):
    with p.open('rb') as stream: digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': p.stat().st_size, 'sha256': digest}
def git(*args, env=None, input=None):
    return subprocess.check_output(['git', *args], env=env, input=input, text=True, timeout=30).strip()
def files(base):
    entries = list(base.rglob('*'))
    assert all(not p.is_symlink() and (p.is_file() or p.is_dir()) for p in entries)
    return sorted(p for p in entries if p.is_file())
def copy_files(base, prefix):
    records = []
    for p in files(base):
        target = ARCHIVE/prefix/p.relative_to(base); target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists(): assert seal(target) == seal(p), p
        else: shutil.copyfile(p, target)
        expected = seal(p); assert seal(target) == expected
        records.append({'origin': str(p.relative_to(ROOT)), 'recovery': str(target.relative_to(ROOT)), **expected})
    return records
def native_seals():
    plan = read(RUNTIME/'execution-plan.json'); prefix = str(NATIVE.relative_to(ROOT))+'/'
    expected = {p:v for p,v in plan['input_seals'].items() if p.startswith(prefix)}
    assert len(expected) == 83
    assert {str(p.relative_to(ROOT)): seal(p) for p in files(NATIVE)} == expected
    return expected
def absent(target):
    pid = int(target['pid']); ticks = str(target['start_ticks']); assert ticks
    try: actual = Path('/proc', str(pid), 'stat').read_text().rsplit(')',1)[1].split()[19]
    except (FileNotFoundError, ProcessLookupError): actual = None
    assert actual != ticks, target
    return {**target, 'recorded_birth_absent': True}

assert len(sys.argv) == 2 and sys.argv[1] in ('archive', 'remove')
assert ROOT.resolve() == Path('/home/astra/.project/Sirius')
assert WORK == ROOT/'out/native-division-first-correction-screen' and not WORK.is_symlink()
assert not git('status','--porcelain') and git('rev-list','--left-right','--count','HEAD...@{u}') == '0\t0'
assert git('rev-list','--left-right','--count','main...origin/main') == '0\t0'
assert git('worktree','list','--porcelain').count('worktree ') == 1
assert git('rev-parse',BRANCH) == BASELINE
assert git('rev-list','--left-right','--count',BRANCH+'...origin/'+BRANCH) == '0\t0'
decision = read(RUNTIME/'root-retention-decision.json')
assert len(decision['retention_checks']) == 44 and len(decision['failed_conditions']) == 33
assert decision['candidate_revision'] == TRIAL and decision['baseline_revision'] == BASELINE
assert read(RUNTIME/'retention-review.json')['accepted'] and read(RUNTIME/'final-ownership-review.json')['accepted']
assert read(RUNTIME/'retention-review.json')['reports'] == read(RUNTIME/'final-ownership-review.json')['reports']
for name, expected in read(RUNTIME/'retention-review.json')['reports'].items(): assert seal(RUNTIME/name) == expected
source_review = read(WORK/'preservation-source-review.json')
assert source_review['accepted'] and source_review['helper_seal'] == seal(Path(__file__))

if sys.argv[1] == 'archive':
    assert git('rev-parse','HEAD') == TRIAL and not ARCHIVE.exists() and not git('tag','--list',TAG)
    for path, expected in read(RUNTIME/'execution-plan.json')['input_seals'].items(): assert seal(ROOT/path) == expected, path
    for path, expected in read(RUNTIME/'matched-native-comparison.json')['complete_input_seals'].items(): assert seal(ROOT/path) == expected, path
    stage = native_seals()
    env = dict(os.environ, GIT_INDEX_FILE=str(WORK/'preservation.index')); blobs = {}
    try:
        git('read-tree', TRIAL, env=env)
        for p in files(WORK):
            relative = p.relative_to(WORK)
            if p.suffix not in ('.py','.ps1','.cs','.slang','.patch') and relative != Path('runtime/native-cases.json'): continue
            if any(part in ('owner-temp','artifact','native-stage-recovery') for part in relative.parts): continue
            if relative.parts[:2] in (('runtime','baseline'),('runtime','source')): continue
            if relative.parts[0] == 'candidate' and p.name != 'retained_triple.slang': continue
            blob = git('hash-object','-w',str(p)); name = PREFIX+str(relative)
            git('update-index','--add','--cacheinfo','100644',blob,name,env=env)
            blobs[str(relative)] = {'git_path': name, 'git_blob': blob}
        tree = git('write-tree',env=env)
        message = (WORK/'preservation-commit-message.txt').read_text()
        assert message.startswith('Preserve rejected native first-correction trial (#44)\n')
        commit = git('commit-tree',tree,'-p',TRIAL,input=message)
        git('tag',TAG,commit)
    finally:
        for p in WORK.glob('preservation.index*'): p.unlink()
    ARCHIVE.mkdir(parents=True)
    records = copy_files(WORK,'files')
    for row in records: row.update(blobs.get(str(Path(row['origin']).relative_to(WORK.relative_to(ROOT))),{}))
    native_records = copy_files(NATIVE,'native-stage')
    assert {r['origin']:{k:r[k] for k in ('bytes','sha256')} for r in native_records} == stage
    manifest = {'source_revision':TRIAL,'restoration_revision':BASELINE,'source_tag':TAG,'source_commit':commit,
        'records':records,'native_stage_records':native_records,'useful_source_paths':len(blobs),
        'historical_versions':{'original_builder':'runtime/baseline/scripts/build-retained-kernels.py',
            'original_triple':'runtime/baseline/src/sirius/kernels/retained_triple.slang',
            'original_header':'runtime/baseline/bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h',
            'candidate_header':'runtime/candidate/retained_kernels.h'},
        'scope':'Full candidate production tree and useful helper/prototype versions plus concise findings in Git. All raw WORK files and all83 stage files locally archived. Original2053 source epoch differs intentionally from runtime ea6 epoch. No remote raw backup/full qualification.'}
    (ARCHIVE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'source_tag':TAG,'source_commit':commit,'work_files':len(records),'native_files':len(native_records),'useful_sources':len(blobs),'bytes':sum(r['bytes'] for r in records+native_records)}))
else:
    assert git('rev-parse','HEAD') == BASELINE and git('branch','--show-current') == BRANCH
    review = read(WORK/'cleanup-ready-review.json'); assert review['accepted'] and review['helper_seal'] == seal(Path(__file__))
    assert read(WORK/'publication-readback.json')['exact']
    restored = read(WORK/'restored-source-build.json')
    assert restored['accepted'] and read(WORK/'restored-products.json')['accepted'] and read(WORK/'restoration-actual-review.json')['accepted']
    bridge = read(WORK/'final-retirement-run.json'); census = read(WORK/'final-native-retirement.json')
    assert bridge['accepted'] and bridge['inputs_unchanged'] and not bridge['errors'] and not bridge['signals']
    assert census['pass'] and not census['matches'] and not census['errors']
    targets = read(WORK/'closing-births.json')
    assert len(census['targets']) == 51 and all(t['recorded_birth_absent'] for t in census['targets'])
    assert [(t['category'],t['kind'],int(t['pid']),t['birth']) for t in census['targets']] == [(t['category'],t['kind'],int(t['pid']),t['birth']) for t in targets['windows_targets']]
    assert len(bridge['saved_linux_births_before']) == len(bridge['saved_linux_births_after']) == 62
    assert all(t['recorded_birth_absent'] for t in bridge['saved_linux_births_before']+bridge['saved_linux_births_after'])
    for path, expected in bridge['input_seals'].items(): assert seal(Path(path)) == expected, path
    assert census['stage'] == subprocess.check_output(['wslpath','-w',str(WORK)],text=True,timeout=10).strip()
    assert census['native_stage'] == subprocess.check_output(['wslpath','-w',str(NATIVE)],text=True,timeout=10).strip()
    manifest = read(ARCHIVE/'manifest.json')
    assert git('ls-remote','origin','refs/tags/'+TAG).split()[0] == manifest['source_commit']
    for row in manifest['records']+manifest['native_stage_records']:
        assert seal(ROOT/row['origin']) == seal(ROOT/row['recovery']) == {k:row[k] for k in ('bytes','sha256')}
        if 'git_blob' in row: assert git('rev-parse',manifest['source_commit']+':'+row['git_path']) == row['git_blob']
    stage = native_seals()
    assert {r['origin']:{k:r[k] for k in ('bytes','sha256')} for r in manifest['native_stage_records']} == stage
    births = [absent(t) for t in targets['linux_targets']]
    births += [absent({'kind':'restored '+r['name'],'pid':r['pid'],'start_ticks':r['start_ticks']}) for r in restored['runs']]
    births.append(absent({'kind':'final retirement observer','pid':bridge['pid'],'start_ticks':bridge['start_ticks']}))
    matches = []; hidden = []; visible = 0
    for p in Path('/proc').iterdir():
        if not p.name.isdigit() or int(p.name) in (os.getpid(),os.getppid()): continue
        try: command = (p/'cmdline').read_bytes().replace(b'\0',b' '); visible += 1
        except (FileNotFoundError,ProcessLookupError,PermissionError) as error:
            hidden.append({'pid':int(p.name),'reason':type(error).__name__}); continue
        if any(v.encode() in command for v in (str(WORK),str(WORK.relative_to(ROOT)),str(NATIVE),str(NATIVE.relative_to(ROOT)))):
            matches.append({'pid':int(p.name),'command':command.decode(errors='replace')})
    assert not matches, matches
    records = copy_files(WORK,'files'); publication = ARCHIVE/'publication'; publication.mkdir()
    (publication/'closing-ledger.json').write_text(json.dumps({'source_tag':TAG,'source_commit':manifest['source_commit'],
        'manifest_seal':seal(ARCHIVE/'manifest.json'),'records':records,'native_stage_records':manifest['native_stage_records'],
        'linux_births':births,'linux_visible_processes':visible,'linux_hidden_records':hidden,
        'scope':'Fresh exact saved births and accessible references; complete local recovery and pushed useful source. No hidden/global/scheduled-consumer absence claim.'},indent=2)+'\n')
    shutil.rmtree(NATIVE); shutil.rmtree(WORK)
    state = {'source_revision':BASELINE,'clean':not git('status','--porcelain'),'tracking':git('rev-list','--left-right','--count','HEAD...@{u}'),
        'main_tracking':git('rev-list','--left-right','--count','main...origin/main'),'worktrees':git('worktree','list'),
        'finished_work_absent':not WORK.exists(),'finished_native_stage_absent':not NATIVE.exists(),
        'out_entries':sorted(p.name for p in (ROOT/'out').iterdir()),'removed_work_files':len(records),'native_stage_files':len(stage),
        'scope':'Rejected trial preserved, accepted source/build products restored. Full science/render/release remains open.'}
    assert state['clean'] and state['tracking'] == state['main_tracking'] == '0\t0'
    (publication/'removal-state.json').write_text(json.dumps(state,indent=2)+'\n'); print(json.dumps(state))
