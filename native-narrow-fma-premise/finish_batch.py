"""Preserve terminal native-narrow FMA premise diagnostics, then remove only exact recoverable scratch."""
from pathlib import Path
import hashlib, json, shutil, subprocess, sys
sys.dont_write_bytecode = True
ROOT = Path('/home/astra/.project/Sirius')
HEAD = 'b81e061fa081cdaed37eca3add950bd254bf296d'
NAMES = ('native-narrow-fma-premise',)
WORKS = [ROOT / 'out' / name for name in NAMES]
TAG = 'evidence/issue-95-native-narrow-fma-premise-b81e061'
ARCHIVE = ROOT / 'attestations/review/b81e061/native-narrow-fma-premise'
PUBLICATION = ROOT / 'attestations/review/b81e061/native-narrow-fma-premise-publication'
CLEANUP = ROOT / 'attestations/cleanup/b81e061/native-narrow-fma-premise-finished'

def run(argv, data=None):
    return subprocess.check_output(argv, cwd=ROOT, input=data)

def query(argv):
    return run(argv).decode().strip()

def seal(path):
    raw = Path(path).read_bytes()
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}

def write(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')

def load(path):
    return json.loads(path.read_text())

def tree_check():
    assert query(['git', 'rev-parse', 'HEAD']) == HEAD
    assert not query(['git', 'status', '--porcelain'])
    assert query(['git', 'rev-list', '--left-right', '--count', 'HEAD...@{upstream}']).split() == ['0', '0']
    assert query(['git', 'worktree', 'list', '--porcelain']).count('worktree ') == 1
    assert all(w == ROOT / 'out' / n and not w.is_symlink() for w, n in zip(WORKS, NAMES))

def inventory(work):
    result = []
    for path in work.rglob('*'):
        assert not path.is_symlink(), path
        if path.is_file(): result.append(path)
        else: assert path.is_dir(), path
    return sorted(result, key=lambda p: str(p.relative_to(work)))

def exact(path, row):
    assert seal(path) == {k: row[k] for k in ('bytes', 'sha256')}, str(path)

def input_checks():
    work=WORKS[0];checked=0;phases=[]
    model=load(work/'model-result.json');compiler=load(work/'compiler-result.json')
    assert model['pass'] and compiler['pass']
    assert model['source_revision']==compiler['source_revision']==HEAD
    assert model['device_workloads']==compiler['device_workloads']==0
    assert model['existing_residual_records']==53524 and model['targeted_records']==3840
    assert model['inputs_unchanged'] and compiler['inputs_unchanged']
    assert load(work/'model-independent-review.json')['accepted']
    assert load(work/'compiled-independent-review.json')['accepted']
    assert load(work/'compiled-inspection.json')['pass']
    assert {p.parent.name for p in work.rglob('owner.json')}=={'model','compiler'}
    for name in ('model','compiler'):
        folder=work/name;owner=load(folder/'owner.json')
        assert owner['status']=='terminal' and owner['passed'] and owner['child_exit']==0
        assert owner['source']=={'revision':HEAD,'status':''}
        assert owner['stop_reason'] is None and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert owner['remaining_owned_processes']==[] and owner['cleanup_errors']==[]
        assert owner['observed_births'] and all(r['absent'] for r in owner['observed_births'])
        before,after=load(folder/'inputs-before.json'),load(folder/'inputs-after.json')
        assert before==after
        action_file=work/('actions.json' if name=='model' else 'compiler-actions.json')
        action=load(action_file)[name]
        assert owner['argv']==action['argv'] and owner['limits']==action['limits']
        for row in before.values():exact(row['resolved_path'],row);checked+=1
        phases.append({'phase':name,'passed':True,'seconds':owner['wall_seconds'],'rss_kib':owner['peak_sampled_rss_kib']})
    before=load(work/'compiler-inputs-before.json')
    assert before==load(work/'compiler-inputs-after.json') and len(before)==81
    for row in model['inputs']:
        exact(ROOT/row['path'],row);checked+=1
    for row in before:exact(row['path'],row);checked+=1
    commands=load(work/'compiler-commands.json');assert len(commands)==20
    for row in commands:
        assert row['status']=='terminal' and row['exit']==0 and row['leader_reaped'] and not row.get('forced_cleanup')
        for key in ('tool','stdout','stderr'):exact(row[key]['path'],row[key]);checked+=1
    for stage in compiler['stages']:
        assert stage['baseline_module_exact'] and stage['original_program_words_exact']
        for artifact in stage['outputs'].values():
            for row in artifact.values():exact(row['path'],row);checked+=1
    return {'references_checked':checked,'terminal_phases':phases,
            'scope':'Finite Fraction model and ordinary native-default SPIR-V compiler preservation only. No native device, complete packed scientific readback, benefit, rendering or qualification result. Baseline official export CI runs independently against pushed Git, not this scratch.'}

def process_scan():
    work=WORKS[0];leaders=[];groups=set();sessions=set()
    for name in ('model','compiler'):
        o=load(work/name/'owner.json')
        assert o['boot_id']==Path('/proc/sys/kernel/random/boot_id').read_text().strip()
        leaders.append((o['controller']['pid'],str(o['controller']['start_ticks'])))
        leaders += [(r['pid'],str(r['saved_start_ticks'])) for r in o['observed_births']]
        sessions.add(o['child_pid'])
    for c in load(work/'compiler-commands.json'):
        p=c['process'];leaders.append((p['pid'],p['birth']));groups.add(p['group'])
        assert p['session'] in sessions
    leaders=sorted(set(leaders));matches=[];raced=0;seen=0
    for d in Path('/proc').iterdir():
        if not d.name.isdecimal():continue
        try:
            s=(d/'stat').read_text().rsplit(')',1)[1].split();pid=int(d.name);birth=s[19];group=int(s[2]);session=int(s[3]);seen+=1
            if (pid,birth) in leaders or group in groups or session in sessions:
                matches.append({'pid':pid,'birth':birth,'group':group,'session':session})
        except (OSError,ValueError,IndexError):raced+=1
    assert not matches,matches
    return {'scope':'Accessible point-in-time recorded PID/birth and inherited compiler groups/owned sessions only; excludes hidden, escaped, arbitrary consumers and scheduled/remote work.','leaders':leaders,'groups':sorted(groups),'sessions':sorted(sessions),'inspected':seen,'raced_or_unreadable':raced,'matches':matches}

def check():
    tree_check(); return {'source_revision':HEAD,'input_checks':input_checks(),'process_scan':process_scan()}

def git_tree(paths):
    children = {}; entries = []
    for name, path in paths:
        if '/' in name:
            first, tail = name.split('/',1); children.setdefault(first,[]).append((tail,path))
        else:
            blob=run(['git','hash-object','-w','--stdin'],path.read_bytes()).decode().strip()
            assert run(['git','cat-file','blob',blob]) == path.read_bytes()
            entries.append((name,f'100644 blob {blob}\t{name}\n'))
    for name, below in children.items(): entries.append((name,f'040000 tree {git_tree(below)}\t{name}\n'))
    return run(['git','mktree'],''.join(line for _,line in sorted(entries)).encode()).decode().strip()

def preserve():
    assert not ARCHIVE.exists()
    post=check();write(WORKS[0]/'terminal-contract.json',post)
    chosen=[]
    keep={'model-result.json','model-independent-review.json','emission.json','compiler-result.json','compiled-inspection.json','compiled-independent-review.json','terminal-contract.json','baseline-ci-dispatch.json'}
    for work in WORKS:
        for p in inventory(work):
            rel=str(p.relative_to(work))
            if p.suffix == '.py' or rel == 'kernels/retained_pair.slang' or rel in keep or (p.parent==work and p.name.endswith('independent-review.json')):
                chosen.append((work.name+'/'+rel,p))
    assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+TAG],cwd=ROOT).returncode == 1
    tree=git_tree(chosen);commit=run(['git','commit-tree',tree],b'Preserve guarded native-default FMA premise findings (#95)\n\nFinite exact model and native compiler preservation only. No device, benefit, rendering or final qualification.\n').decode().strip()
    run(['git','tag',TAG,commit]);run(['git','push','origin','refs/tags/'+TAG])
    assert query(['git','ls-remote','--tags','origin','refs/tags/'+TAG]).split() == [commit,'refs/tags/'+TAG]
    source={'tag':TAG,'commit':commit,'remote_exact':True,'diagnostic_only':True,'files':[{'name':n,**seal(p)} for n,p in chosen]};write(WORKS[0]/'source-preservation.json',source)
    ARCHIVE.mkdir(parents=True);entries=[]
    for work in WORKS:
        for p in inventory(work):
            rel=Path(work.name)/p.relative_to(work);dest=ARCHIVE/'diagnostic'/rel;dest.parent.mkdir(parents=True,exist_ok=True);s=seal(p);shutil.copyfile(p,dest);exact(dest,s)
            entries.append({'path':str(dest.relative_to(ARCHIVE)),'scratch_path':str(Path('out')/rel),**s})
    manifest={'source_revision':HEAD,'source_tag':TAG,'source_snapshot':commit,'files':entries,'payload_count':len(entries),'payload_bytes':sum(r['bytes'] for r in entries),'scope':post['input_checks']['scope']};write(ARCHIVE/'manifest.json',manifest)
    print(json.dumps({'archive':str(ARCHIVE.relative_to(ROOT)),'manifest':seal(ARCHIVE/'manifest.json'),'source_snapshot':commit,'source_tag':TAG,'source_files':len(chosen),'payload_count':len(entries),'payload_bytes':manifest['payload_bytes']}))

def cleanup():
    tree_check();assert not CLEANUP.exists()
    review=load(PUBLICATION/'cleanup-independent-review.json');assert review['accepted']
    assert review['manifest']==seal(ARCHIVE/'manifest.json') and review['reviewed_finisher']==seal(Path(__file__))
    manifest=load(ARCHIVE/'manifest.json');source=load(WORKS[0]/'source-preservation.json')
    assert source['commit']==manifest['source_snapshot'] and source['tag']==TAG
    assert query(['git','ls-remote','--tags','origin','refs/tags/'+TAG]).split()==[source['commit'],'refs/tags/'+TAG]
    actual=sorted(str(p.relative_to(ROOT)) for w in WORKS for p in inventory(w));assert actual==sorted(r['scratch_path'] for r in manifest['files'])
    for r in manifest['files']: exact(ROOT/r['scratch_path'],r);exact(ARCHIVE/r['path'],r)
    for r in source['files']:
        raw=run(['git','show',source['commit']+':'+r['name']]);assert {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}=={k:r[k] for k in ('bytes','sha256')}
    terminal=check()
    protected=[ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h',*sorted((ROOT/'bin/linux-gcc/src/sirius/backend/retained').glob('*.spv')),ROOT/'bin/linux-gcc/tests/backend/sirius_backend_tests']
    before={str(p.relative_to(ROOT)):seal(p) for p in protected};assert len(before)==35
    for work in WORKS: shutil.rmtree(work)
    assert not any(w.exists() for w in WORKS)
    assert before=={str(p.relative_to(ROOT)):seal(p) for p in protected};tree_check()
    CLEANUP.mkdir(parents=True);receipt={'source_revision':HEAD,'removed_only':[str(w.relative_to(ROOT)) for w in WORKS],'removed_regular_files':len(manifest['files']),'removed_bytes':manifest['payload_bytes'],'exact_recoveries':manifest['files'],'manifest':seal(ARCHIVE/'manifest.json'),'source_tag':TAG,'source_snapshot':source['commit'],'terminal_contract':terminal,'protected_build_payloads_unchanged':before,'tracked_tree_clean_upstream_zero_zero_single_worktree':True};write(CLEANUP/'cleanup.json',receipt)
    print(json.dumps({'cleanup':str((CLEANUP/'cleanup.json').relative_to(ROOT)),'removed_regular_files':len(manifest['files']),'removed_bytes':manifest['payload_bytes'],'receipt':seal(CLEANUP/'cleanup.json')}))

if __name__=='__main__':
    assert len(sys.argv)==2 and sys.argv[1] in ('--check','--preserve','--cleanup')
    if sys.argv[1]=='--check': print(json.dumps(check(),indent=2))
    elif sys.argv[1]=='--preserve': preserve()
    else: cleanup()
