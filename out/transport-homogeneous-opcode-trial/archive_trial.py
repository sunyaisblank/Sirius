"""Preserve and remove only the finished, rejected issue85 workspace/stage."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, os, shutil, subprocess, sys

ROOT=Path.cwd()
WORK=Path(__file__).resolve().parent
ARCHIVE=ROOT/'attestations/native-vulkan/bd82e52/homogeneous-opcode-dispatch-rejected'
PROOF=ROOT/'attestations/cleanup/2026-10-10/homogeneous-opcode-dispatch-final-proofs'
RECEIPT=ROOT/'attestations/cleanup/2026-10-10/homogeneous-opcode-dispatch-finished.json'
HEAD='90e1de093d419dfa1e78f15544540487bca0dbad'
STAGE=ROOT/'bin/windows-msvc/native-bd82e52'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def doc(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def save(p,d): p.write_text(json.dumps(d,indent=2)+'\n')
def git(*args): return subprocess.check_output(['git',*args],text=True).strip()
def verify_git():
    assert git('rev-parse','HEAD')==HEAD and not git('status','--porcelain')
    assert git('rev-list','--left-right','--count','HEAD...@{upstream}')=='0\t0'
    assert git('rev-parse','HEAD^{tree}')=='a336cf8544a127a350da7fac0e5a5aed75e5503b'
    assert git('ls-remote','origin','refs/heads/development/retained-dense-step-control-2026-10-06').split()[0]==HEAD
def verify_recovery():
    data=doc(WORK/'whole-native-stage-and-zip-recovery.json')
    for record in data['stage_maps'].values():
        paths={str(p.relative_to(ROOT)) for p in (ROOT/record['stage']).rglob('*') if p.is_file()}
        assert paths=={r['stage_path'] for r in record['payloads']}
        for r in record['payloads']:
            for rel in [r['stage_path'],*r['retained_export_paths']]:
                p=ROOT/rel
                assert p.is_file() and not p.is_symlink() and p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],rel
    for rel,r in data['zip_payload_maps'].items():
        p=ROOT/rel
        assert p.stat().st_size==r['container_bytes'] and sha(p)==r['container_sha256']
        for q in r['payloads']:
            p=ROOT/q['retained_path']
            assert p.stat().st_size==q['bytes'] and sha(p)==q['sha256']
    return data
def verify_source_preservation():
    d=doc(WORK/'complete-tool-source-git-preservation.json')
    assert d['remote_tag_and_commit_verified'] and d['canonical_HEAD_index_worktree_unchanged']
    commit=d['source_commit']
    assert git('ls-remote','origin','refs/tags/'+d['tag']+'^{}').split()[0]==commit
    hashes=set()
    for r in d['files']:
        p=ROOT/r['path']
        assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256']
        raw=subprocess.check_output(['git','show',commit+':'+r['path']])
        assert raw==p.read_bytes()
        hashes.add(r['sha256'])
    for p in WORK.rglob('*'):
        if p.is_file() and p.suffix in ('.py','.ps1'):
            assert sha(p) in hashes,str(p)
    return d
def fields(pid): return Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()
def linux_final(data):
    births=[]
    for r in data['linux_births']:
        try: actual=fields(r['pid'])[19]
        except (FileNotFoundError,ProcessLookupError): actual=None
        assert actual!=r['start_ticks'],r
        births.append({'pid':r['pid'],'saved_start_ticks':r['start_ticks'],'current_start_ticks':actual,'recorded_birth_absent':True})
    ancestors=[];pid=os.getpid()
    while pid>0:
        try:f=fields(pid)
        except (FileNotFoundError,ProcessLookupError):break
        ancestors.append({'pid':pid,'start_ticks':f[19]});pid=int(f[1])
    excluded={(r['pid'],r['start_ticks']) for r in ancestors}
    matches=[];unreadable={};count=0
    for folder in Path('/proc').iterdir():
        if not folder.name.isdecimal():continue
        pid=int(folder.name)
        try:f=fields(pid)
        except (FileNotFoundError,ProcessLookupError):continue
        except OSError:unreadable['stat']=unreadable.get('stat',0)+1;continue
        if (pid,f[19]) in excluded:continue
        count+=1;hits=[]
        for name in ['cmdline','exe','cwd']:
            try:
                value=(folder/name).read_bytes().replace(b'\0',b' ').decode(errors='replace') if name=='cmdline' else os.readlink(folder/name)
                if any(word in value for word in data['linux_targets']):hits.append(name)
            except (FileNotFoundError,ProcessLookupError):pass
            except OSError:unreadable[name]=unreadable.get(name,0)+1
        try:
            for fd in (folder/'fd').iterdir():
                try:
                    if any(word in os.readlink(fd) for word in data['linux_targets']):hits.append('fd:'+fd.name)
                except (FileNotFoundError,ProcessLookupError):pass
                except OSError:unreadable['fd_target']=unreadable.get('fd_target',0)+1
        except (FileNotFoundError,ProcessLookupError):pass
        except OSError:unreadable['fd_directory']=unreadable.get('fd_directory',0)+1
        if hits:matches.append({'pid':pid,'start_ticks':f[19],'fields':hits})
    assert not matches,matches
    return {'saved_births':births,'recorded_births_verified':len(births),'processes':count,'matches':matches,'unreadable_fields':unreadable,'excluded_exact_inspector_ancestors':ancestors,'scope':'Exact recorded Linux births and accessible cmdline/exe/cwd/fd removal-target references; unreadable/global consumers not excluded.'}

mode,=sys.argv[1:];assert mode in ('archive','remove')
verify_git();recovery=verify_recovery();source=verify_source_preservation()
assert doc(WORK/'restoration-completion.json')['pass']
assert doc(WORK/'independent-native-retention-review.json')['failed_checks']==8
assert not doc(WORK/'native-retention-comparison.json')['retention_gate_pass']
if mode=='archive':
    assert not ARCHIVE.exists();ARCHIVE.mkdir(parents=True)
    exclusions={str(Path(rel).relative_to(WORK.relative_to(ROOT))) for rel in recovery['zip_payload_maps']}
    records=[]
    for p in sorted(WORK.rglob('*')):
        assert not p.is_symlink()
        if not p.is_file():continue
        rel=str(p.relative_to(WORK))
        if rel in exclusions:continue
        r={'path':rel,'bytes':p.stat().st_size,'sha256':sha(p)}
        target=ARCHIVE/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
        assert target.stat().st_size==r['bytes'] and sha(target)==r['sha256']
        records.append(r)
    manifest={'candidate_revision':'bd82e523227a1b7223f43984af10291711925e44','baseline_revision':'eef7baf23f2e8c9906a467a6eef9ba7b6317d5ee','restoration_revision':HEAD,'files':records,'files_count':len(records),'bytes':sum(r['bytes'] for r in records),'excluded_zip_containers':recovery['zip_payload_maps'],'tool_source_commit':source['source_commit'],'tool_source_tag':source['tag'],'disposition':'rejected;8/46 gate failures;accepted production restored;no retries','scope':'Whole local ignored trial evidence/source payload archive. All ZIP members remain in verified official bundles; compressed containers are excluded and not claimed byte-recreated. No remote raw-evidence backup/fullqualification claim.'}
    save(ARCHIVE/'manifest.json',manifest)
    print(json.dumps({'archive_files':len(records),'archive_bytes':manifest['bytes'],'manifest_sha256':sha(ARCHIVE/'manifest.json')}))
else:
    assert not RECEIPT.exists()
    manifest=doc(ARCHIVE/'manifest.json')
    expected={r['path'] for r in manifest['files']}
    exclusions=set(str(Path(rel).relative_to(WORK.relative_to(ROOT))) for rel in recovery['zip_payload_maps'])
    assert {str(p.relative_to(WORK)) for p in WORK.rglob('*') if p.is_file()}==expected|exclusions
    assert {str(p.relative_to(ARCHIVE)) for p in ARCHIVE.rglob('*') if p.is_file()}==expected|{'manifest.json'}
    for r in manifest['files']:
        for folder in [WORK,ARCHIVE]:
            p=folder/r['path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256']
    final=doc(PROOF/'cleanup-trial-accepted.json');bridge=doc(PROOF/'cleanup-trial-bridge.json');windows=doc(PROOF/'cleanup-trial.json')
    assert final['pass'] and bridge['accepted'] and bridge['linux_child_birth_absent'] and not bridge['errors'] and not bridge['outer_stop']
    assert windows['pass'] and windows['known_births_verified']==46 and not windows['matches'] and not windows['accessible_compilers']
    assert sha(PROOF/'cleanup-trial.json')==final['terminal_receipt_sha256']
    assert sha(PROOF/'cleanup-trial-bridge.json')==final['terminal_bridge_sha256']
    inputs=doc(WORK/'trial-final-cleanup-inputs.json');linux=linux_final(inputs)
    save(PROOF/'linux-final-consumers.json',linux)
    for name in ['archive_trial.py','run_trial_final_audit.py','trial-final-cleanup-audit.ps1','trial-final-cleanup-inputs.json','trial-final-cleanup-owner-review.json','whole-native-stage-and-zip-recovery.json']:
        shutil.copyfile(WORK/name,PROOF/name);assert sha(WORK/name)==sha(PROOF/name)
    # All source/payload, official recovery and publisher joins remain exact immediately before root-only removal.
    verify_git();verify_recovery();verify_source_preservation()
    for r in manifest['files']:
        for folder in [WORK,ARCHIVE]:
            p=folder/r['path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256']
    publication=doc(WORK/'candidate-rejection-publication.json')
    assert publication['readback']['body']==(WORK/'candidate-rejection-publication.txt').read_text()
    result={'completed_utc':datetime.now(timezone.utc).isoformat(),'source_revision':HEAD,'branch':'development/retained-dense-step-control-2026-10-06','tracking':'0/0','committed_and_pushed':True,'archive':str(ARCHIVE.relative_to(ROOT)),'archive_payload_files':manifest['files_count'],'archive_payload_bytes':manifest['bytes'],'manifest_sha256':sha(ARCHIVE/'manifest.json'),'whole_payloads_reverified':True,'tool_source_commit':source['source_commit'],'tool_source_tag':source['tag'],'rejection_comment':publication['readback']['html_url'],'windows':windows,'linux':linux,'rejected_stage':recovery['stage_maps']['candidate'],'removed_targets':[str(WORK.relative_to(ROOT)),str(STAGE.relative_to(ROOT))],'finished_scratch_removed':False,'rejected_stage_removed':False,'official_bundles_retained':4,'accepted_baseline_stage_retained':True,'raw_evidence_remote_backup_claimed':False,'qualification_claimed':False}
    save(RECEIPT,result)
    shutil.rmtree(STAGE);assert not STAGE.exists()
    result['rejected_stage_removed']=True;save(RECEIPT,result)
    shutil.rmtree(WORK);assert not WORK.exists()
    result['finished_scratch_removed']=True;save(RECEIPT,result)
    print(json.dumps({'removed_targets':result['removed_targets'],'archive_files':manifest['files_count'],'archive_bytes':manifest['bytes'],'linux_births_verified':linux['recorded_births_verified'],'windows_births_verified':windows['known_births_verified']}))
