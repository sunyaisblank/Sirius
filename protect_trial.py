"""Preserve a finished rejected trial before restoring any volatile products."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, os, shutil, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p): return {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
def key(v): return v['bytes'], v['sha256']
def expected(v): return {k:v[k] for k in ('bytes','sha256')}
def git(*args, **kw): return subprocess.check_output(['git', *args], **kw)
source = read(WORK/'source.json')
head, baseline = source['source_revision'], source['baseline_revision']
assert git('rev-parse','HEAD',text=True).strip()==head and not git('status','--porcelain')
assert git('rev-parse','HEAD')==git('rev-parse','@{u}')
disposition = read(WORK/'native-disposition-reviewed.json')
assert disposition['accepted'] and disposition['batch_terminal'] and disposition['reject_candidate']
assert disposition['source_revision']==head and disposition['baseline_revision']==baseline
comparison = read(WORK/'native-retention-comparison.json')
assert not comparison['retention_gate_pass'] and len(comparison['checks'])==46 and len(comparison['failed_checks'])==25
assert disposition['comparison_sha256']==seal(WORK/'native-retention-comparison.json')['sha256']
controls=read(WORK/'native-controls-analysis.json')
assert controls['pass'] and controls['total_cases']==10 and len(controls['results'])==9 and len(controls['additive_controls'])==1
for mode,calls,runs in (('controls',20,10),('matched',8,4)):
    sequence=read(WORK/('native-'+mode+'-sequence.json'))
    assert sequence['completed'] and not sequence.get('failure') and len(sequence['results'])==runs and len(sequence['owned_calls'])==calls
    assert all(c['returncode']==0 and c['ownership_verified'] and not c['cleanup_errors'] and not c['controller_signals'] and not c['relay_errors'] for c in sequence['owned_calls'])
for origin,v in read(WORK/'execution-bindings.json')['files'].items():
    assert seal(ROOT/origin)==expected(v),origin
assert seal(WORK/'compare_native.py')['sha256']=='7e13981d4fb3e566c5e5ae404946a899530096d0d3c6f7b09c664fc5f8a5f686'
DEST = ROOT/'attestations/native-vulkan'/head[:7]/'guarded-native-default-fma-rejected'
assert not DEST.exists()
TAG = 'evidence/issue-95-native-narrow-fma-trial-'+head[:7]
assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+TAG]).returncode==1

# Source-only Git tree: used helpers, exact historical helper versions, fixed
# controls and concise reviewed findings. Never replace the checkout's index.
names = {p for p in WORK.glob('*.py')} | {p for p in WORK.glob('*.ps1')}
names |= set((WORK/'historical-inputs').rglob('*'))
names |= set(WORK.glob('*actions.json'))
names |= {WORK/n for n in ('source.json','candidate.patch','retention-gate-frozen.json',
    'native-plan.json','expected-readback-layout.json','independent-source-stream-review.json',
    'independent-build-admission-review.json','baseline-native-independent-review.json',
    'candidate-native-apparatus-review.json','native-controls-independent-review.json',
    'native-disposition-reviewed.json','protection-source-review.json')}
names = sorted(p for p in names if p.is_file())
index = WORK/'source-tools.index'
assert not index.exists()
env = {**os.environ,'GIT_INDEX_FILE':str(index)}
subprocess.run(['git','read-tree','--empty'],env=env,check=True)
tools=[]
for p in names:
    assert not p.is_symlink()
    name=str(p.relative_to(WORK));blob=git('hash-object','-w','--no-filters',str(p),text=True).strip()
    subprocess.run(['git','update-index','--add','--cacheinfo','100644',blob,name],env=env,check=True)
    tools.append({'path':name,'origin':str(p.relative_to(ROOT)),'blob':blob,**seal(p)})
tree=git('write-tree',env=env,text=True).strip()
commit=git('commit-tree',tree,input=(f'Preserve rejected guarded native-default FMA trial (#95)\n\nCandidate {head}; baseline {baseline}. Ten native controls passed; original46 gate failed 25. No rendering or full qualification.\n').encode()).decode().strip()
subprocess.run(['git','tag',TAG,commit],check=True)
subprocess.run(['git','push','origin','refs/tags/'+TAG],check=True,timeout=90)
assert git('ls-remote','origin','refs/tags/'+TAG,text=True).split()[0]==commit
index.unlink()
preserved={'tag':TAG,'commit':commit,'source_revision':head,'baseline_revision':baseline,
    'pushed':True,'remote_commit_verified':True,'files':tools}
(WORK/'tools-preservation.json').write_text(json.dumps(preserved,indent=2)+'\n')

recoveries={}
def add(value,recovery): recoveries.setdefault(key(value),recovery)
# Index immutable production blobs with their actual bytes, including old model
# identities. A content match alone does not claim the original path was current.
blobs={}
for revision in (head,baseline):
    for row in git('ls-tree','-r','-z',revision).split(b'\0'):
        if not row: continue
        info,path=row.split(b'\t',1);mode,kind,oid=info.decode().split();assert kind=='blob'
        blobs.setdefault(oid,{'kind':'immutable-Git-blob','revision':revision,'blob':oid,'path':path.decode()})
for t in tools: blobs[t['blob']]={'kind':'immutable-Git-blob','revision':commit,'blob':t['blob'],'path':t['path']}
oids=list(blobs);data=git('cat-file','--batch',input=('\n'.join(oids)+'\n').encode());offset=0
for oid in oids:
    end=data.index(b'\n',offset);header=data[offset:end].split();assert header[0].decode()==oid and header[1]==b'blob'
    n=int(header[2]);body=data[end+1:end+1+n];assert data[end+1+n:end+2+n]==b'\n'
    add({'bytes':n,'sha256':hashlib.sha256(body).hexdigest()},blobs[oid]);offset=end+2+n
assert offset==len(data)
del data
for revision in (baseline,head):
    for p in (ROOT/'attestations/native-build'/revision[:7]).rglob('*'):
        assert not p.is_symlink()
        if p.is_file(): add(seal(p),{'kind':'retained-native-export','path':str(p.relative_to(ROOT))})
for relative in ('attestations/native-vulkan/95c36cb/original-fp64-observation',
                 'attestations/review/b81e061/native-narrow-fma-premise',
                 'attestations/review/b81e061/compact-transport-feasibility',
                 'attestations/review/b81e061/literal-transport-compiler-feasibility'):
    archive=ROOT/relative;manifest=read(archive/'manifest.json')
    for v in manifest.get('payloads',[]):
        if 'bytes' in v and 'sha256' in v: add(v,{'kind':'retained-prior-payload','path':str((archive/v['path']).relative_to(ROOT))})
    for v in manifest.get('executed_input_ledger',[]):
        if 'recovery' in v: add(v,v['recovery'])

payloads={}
def resolve(path,value):
    k=key(value)
    if k not in recoveries:
        assert path.is_file() and seal(path)==expected(value),str(path)
        target='payloads/'+value['sha256'];dst=DEST/target;dst.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,dst);assert seal(dst)==expected(value)
        payloads[target]={'path':target,**expected(value)}
        add(value,{'kind':'local-archive-payload','path':str(dst.relative_to(ROOT))})
    r=recoveries[k]
    if r['kind']=='immutable-Git-blob':
        b=git('cat-file','blob',r['blob']);assert {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}==expected(value)
    else: assert seal(ROOT/r['path'])==expected(value),r
    return r

references={}
for p in [*WORK.glob('*/inputs-before.json'),*WORK.glob('*/inputs-after.json'),
          WORK/'baseline-native-bindings.json',WORK/'execution-bindings.json']:
    if p.is_relative_to(WORK/'protection'): continue
    d=read(p);d=d.get('files',d)
    for origin,v in d.items():
        if isinstance(v,dict) and isinstance(v.get('bytes'),int) and isinstance(v.get('sha256'),str):
            entry=references.setdefault((origin,*key(v)),{'origin':origin,**expected(v),'ledgers':[]})
            entry['ledgers'].append(str(p.relative_to(ROOT)))
ledger=[]
for (origin,n,digest),v in sorted(references.items()):
    p=ROOT/origin;current=p.is_file() and seal(p)==expected(v)
    # Historical versions resolve by exact pushed tool bytes, never by pretending
    # a changed private helper stayed current across completed owners.
    ledger.append({**v,'origin_current_at_preservation':current,'recovery':resolve(p,v)})
historical={(v['origin'],v['bytes'],v['sha256']) for v in ledger if not v['origin_current_at_preservation']}
assert historical=={
    ('out/native-narrow-fma-trial/reconstruct_native.py',5668,'739db8ec36bd4bc09c48352040c272ac77b9dff8673708bcc08dd565babadf44'),
    ('out/native-narrow-fma-trial/reconstruct_native.py',5186,'70cc4c66c482ad1df97154da85cd756e9edfbda4cba4fb81b9e98a2580ab7f0d'),
    ('out/native-narrow-fma-trial/download_exports.py',3420,'55cf636f640f04f9713508286df961b0aceddcab61c4ac262ad576c150f893bf')}

diagnostics=[];additional=[]
for p in sorted(WORK.rglob('*')):
    assert not p.is_symlink()
    if p.is_relative_to(WORK/'protection'): continue
    if p.is_file():
        v=seal(p);diagnostics.append({'origin':str(p.relative_to(ROOT)),**v,'recovery':resolve(p,v)})
paths=[ROOT/'bin/linux-gcc/src/sirius/backend/retained',ROOT/'bin/linux-gcc/generated',
       ROOT/'bin/linux-gcc/CMakeCache.txt',ROOT/'bin/linux-gcc/build.ninja']
for platform in ('windows-msvc','macos-clang'):
    paths += [ROOT/'bin'/platform/('native-'+r[:7]) for r in (baseline,head)]
for directory in paths:
    assert directory.exists() and not directory.is_symlink()
    for p in sorted(directory.rglob('*')) if directory.is_dir() else [directory]:
        assert not p.is_symlink()
        if p.is_file():
            v=seal(p);additional.append({'origin':str(p.relative_to(ROOT)),**v,'recovery':resolve(p,v)})
manifest={'source_revision':head,'baseline_revision':baseline,'disposition':disposition,
    'tools':preserved,'executed_input_ledger':ledger,'diagnostics':diagnostics,
    'volatile_product_and_native_stage_files':additional,'payloads':list(payloads.values()),
    'payload_count':len(payloads),'payload_bytes':sum(v['bytes'] for v in payloads.values()),
    'scope':'Exact finished rejected native trial inputs and raw outputs, with three historical helper-version differences resolved by exact pushed bytes. Four original API-digest ZIPs and all four retained official native exports remain recoverable. Finite native numerical/performance evidence; sampled RSS, uncalibrated Windows/WSL clocks and bounded accessible ownership visibility. No frame/full-science/release qualification or cleanup authorization.'}
DEST.mkdir(parents=True,exist_ok=True)
(DEST/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
assert git('rev-parse','HEAD',text=True).strip()==head and not git('status','--porcelain')
receipt={'archive':str(DEST.relative_to(ROOT)),'manifest':seal(DEST/'manifest.json'),
    'payload_count':len(payloads),'payload_bytes':manifest['payload_bytes'],
    'executed_origin_version_pairs':len(ledger),'historical_versions':sum(not v['origin_current_at_preservation'] for v in ledger),
    'diagnostic_files':len(diagnostics),'additional_product_stage_files':len(additional),
    'source_tag':TAG,'source_commit':commit,'all_recovery_bytes_verified':True,'cleanup_authorized':False}
(WORK/'trial-protection.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
