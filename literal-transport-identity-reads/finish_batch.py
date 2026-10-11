"""Preserve four terminal diagnostics, then remove only their exact recoverable scratch."""
from pathlib import Path
import hashlib, json, shutil, subprocess, sys
sys.dont_write_bytecode = True
ROOT = Path('/home/astra/.project/Sirius')
HEAD = 'b81e061fa081cdaed37eca3add950bd254bf296d'
NAMES = ('literal-transport-identity-reads', 'literal-transport-table-storage',
         'literal-transport-recognizer-isolation', 'literal-transport-recognizer-marker')
WORKS = [ROOT / 'out' / name for name in NAMES]
TAG = 'evidence/issue-44-transport-compiler-feasibility-b81e061'
ARCHIVE = ROOT / 'attestations/review/b81e061/literal-transport-compiler-feasibility'
PUBLICATION = ROOT / 'attestations/review/b81e061/literal-transport-feasibility-publication'
CLEANUP = ROOT / 'attestations/cleanup/b81e061/literal-transport-feasibility-finished'

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

def compiler_refs(work):
    before, result = load(work/'compile-before.json'), load(work/'compile-result.json')
    assert before['revision'] == HEAD and result['pass'] and result['source_unchanged'] and result['stage_inputs_unchanged']
    refs = [before[k] for k in ('source','emitter','runner','program_metadata','emission_metadata','schedule','header','builder','baseline','baseline_copy')]
    refs += before['stage_inputs'] + list(before['tools'].values())
    emission = load(work/'emission.json')
    refs += [{**r, 'path': str(ROOT/r['path'])} for r in emission['inputs']+[emission['header']]]
    commands = load(work/'commands.json')
    assert len(commands) == 7
    for record in commands:
        assert record['disposition'] == 'completed' and record['returncode'] == 0
        c = record['cleanup']
        assert c['known_leader_reaped'] and not c['errors'] and not c['forced_cleanup']
        assert c['terminal_group_scan'] is not None and not c['terminal_group_scan']['members']
        refs += [record['tool'], record['stdout_seal'], record['stderr_seal']]
    refs += [result['module'], result['assembly']]
    assert load(work/'compiled-inspection.json')['pass']
    return refs

def input_checks():
    checked, versions, phases = 0, [], []
    expected_phases = [set(), {'prepare-hosts','candidate-timestamps','join-original-producer'}, {'prepare-hosts','candidate-timestamps'}, {'prepare-hosts','candidate-timestamps'}]
    for work, expected in zip(WORKS, expected_phases):
        discovered = {str(p.parent.relative_to(work)) for p in work.rglob('owner.json')}
        assert discovered == expected, (work.name, discovered, expected)
        assert (work/'host-preparation.json').exists() is bool(expected)
        for name in expected:
            assert all((work/name/n).is_file() for n in ('owner.json','inputs-before.json','inputs-after.json','stdout.log','stderr.log','samples.jsonl')), (work.name,name)
        refs = compiler_refs(work)
        if (work/'host-preparation.json').exists():
            h = load(work/'host-preparation.json')
            assert h['pass'] and h['device_workloads'] == 0 and h['all_inputs_unchanged'] and h['input_count'] == 949
            before, after = load(work/'host-inputs-before.json'), load(work/'host-inputs-after.json')
            assert before == after
            refs += list(before.values())
            for a in h['artifacts'].values(): refs += [a['executable'], a['factory']]
            assert all(c['exit'] == 0 for c in load(work/'host-commands.json'))
            assert load(work/'host-payload-check.json')['pass']
        for folder in sorted(work.iterdir()):
            if not folder.is_dir() or not (folder/'owner.json').exists(): continue
            owner = load(folder/'owner.json')
            assert owner['status'] == 'terminal' and owner['source'] == {'revision':HEAD,'status':''}
            assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
            assert owner['remaining_owned_processes'] == [] and owner['cleanup_errors'] == []
            rejected = work.name == NAMES[1] and folder.name == 'candidate-timestamps'
            assert owner['passed'] is (not rejected)
            if rejected:
                assert owner['child_exit'] == -15 and owner['stop_reason'] == 'diagnostic_sampled_rss_guard'
                assert not (folder/'gtest.xml').exists()
            else:
                assert owner['child_exit'] == 0 and owner['stop_reason'] is None
                if folder.name == 'candidate-timestamps':
                    from xml.etree import ElementTree as ET
                    x=ET.parse(folder/'gtest.xml').getroot()
                    assert x.tag=='testsuites' and x.attrib['tests']=='1'
                    assert all(x.attrib[k]=='0' for k in ('failures','errors','disabled'))
                    cases=x.findall('./testsuite/testcase');assert len(cases)==1
                    assert cases[0].attrib['name']=='DeviceTimestampsPreserveOriginalIntervalResults' and not cases[0].findall('skipped')
            before, after = load(folder/'inputs-before.json'), load(folder/'inputs-after.json')
            assert before == after
            version = work/'action-versions'/(folder.name+'.json')
            action_key = str((work/'actions.json').relative_to(ROOT))
            exact(version, before[action_key]); versions.append({'phase':str(folder.relative_to(ROOT)), 'version':str(version.relative_to(ROOT)), **seal(version)})
            for name, row in before.items():
                path = version if name == action_key else Path(row['resolved_path'])
                exact(path, row); checked += 1
            phases.append({'phase':str(folder.relative_to(ROOT)), 'passed':owner['passed'], 'seconds':owner['wall_seconds'], 'rss_kib':owner['peak_sampled_rss_kib'], 'stop_reason':owner['stop_reason']})
        for row in refs:
            exact(row['path'], row); checked += 1
    producer = load(WORKS[1]/'original-producer-join.json')
    assert producer['pass'] and producer['exit'] == 0 and producer['inputs_before'] == producer['inputs_after']
    for name, row in producer['inputs_before'].items(): exact(name, row); checked += 1
    exact(producer['relinked_product']['path'], producer['relinked_product'])
    assert seal(producer['relinked_product']['path']) == {k:producer['accepted_product'][k] for k in ('bytes','sha256')}
    for work in WORKS[2:]:
        r = load(work/'numerical-control-readback.json')
        assert r['pass'] and r['full_readback_digests_exact'] == 32 and r['ordered_available_markers'] == 51
        assert r['historical_baseline_producer_revision'] == '7f028007a4935ffdefe5ee9c886b138299f50dd8'
    return {'references_checked':checked, 'exact_historical_action_versions':versions, 'terminal_phases':phases,
            'scope':'Exact recorded inputs and positive compiler/build contracts; bounded full-candidate RSS rejection and two generic numerical controls. No full literal numerical, performance, rendering, adoption or final qualification.'}

def process_scan():
    leaders, groups, sessions = [], set(), set()
    for work in WORKS:
        o = load(work/'compile-before.json')['owner']; leaders.append((o['pid'], str(o['birth'])))
        for c in load(work/'commands.json'):
            leaders.append((c['pid'], str(c['birth']))); groups.add(c['pid'])
        for p in work.rglob('owner.json'):
            o = load(p); assert o['boot_id'] == Path('/proc/sys/kernel/random/boot_id').read_text().strip()
            leaders.append((o['controller']['pid'], str(o['controller']['start_ticks'])))
            assert all(r['absent'] is True for r in o['observed_births'])
            leaders += [(r['pid'], str(r['saved_start_ticks'])) for r in o['observed_births']]
            sessions.add(o['child_pid'])
    matches, raced, seen = [], 0, 0
    for d in Path('/proc').iterdir():
        if not d.name.isdecimal(): continue
        try:
            s=(d/'stat').read_text().rsplit(')',1)[1].split(); pid=int(d.name); birth=s[19]; group=int(s[2]); session=int(s[3]); seen+=1
            if (pid,birth) in leaders or group in groups or session in sessions:
                matches.append({'pid':pid,'birth':birth,'group':group,'session':session})
        except (OSError, ValueError, IndexError): raced+=1
    assert not matches, matches
    return {'scope':'Accessible point-in-time /proc recorded PID/birth, compiler groups and owned host/runtime sessions; excludes hidden, escaped and scheduled work', 'leaders':leaders,'groups':sorted(groups),'sessions':sorted(sessions),'inspected':seen,'raced_or_unreadable':raced,'matches':matches}

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
    keep={'batch-disposition.json','emission.json','programs.json','schedule.json','compile-result.json','compiled-inspection.json','source-frame-check.json','identity-source-check.json','terminal-contract.json','numerical-control-readback.json','feasibility-disposition.json','original-producer-join.json','marker-proposal.json','isolation-proposal.json','proposal.json','recognizer-frame-check.json'}
    for work in WORKS:
        for p in inventory(work):
            rel=str(p.relative_to(work))
            if p.suffix == '.py' or rel in keep or (p.parent==work and p.name.endswith('independent-review.json')):
                chosen.append((work.name+'/'+rel,p))
    assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+TAG],cwd=ROOT).returncode == 1
    tree=git_tree(chosen);commit=run(['git','commit-tree',tree],b'Preserve bounded Transport compiler and recognizer feasibility findings (#44)\n\nUnadopted diagnostic source snapshot; full literal candidate exceeded RSS bound. No rendering or final qualification.\n').decode().strip()
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
