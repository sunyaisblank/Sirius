"""Preserve terminal finite evidence and authentic products with verified blob reuse."""
from pathlib import Path
import hashlib, json, os, stat, subprocess, zipfile
R=Path.cwd(); W=Path(__file__).resolve().parent
SOURCE=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
D=R/'attestations/qualification'/SOURCE[:8]/'lavapipe-owner-machine-codegen-finite-outcome'
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
def write(p,x):
    assert not p.exists(),p
    p.write_text(json.dumps(x,indent=2)+'\n')
assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert not subprocess.check_output(['git','diff','0109ab462edade524667d20303507aa09adcbaf6','--','src','scripts','tests','.github','CMakeLists.txt','CMakePresets.json']).strip()
assert not D.exists()
owners={}; products=set(); sealed_build={}
for p in sorted(W.rglob('owner.json')):
    o=doc(p)
    assert o['status']=='terminal' and not o['cleanup_errors'] and not o['remaining_owned_processes']
    assert o['source_unchanged'] and o['whole_inputs_unchanged']
    assert o['observed_births'] and all(x['absent'] for x in o['observed_births'])
    for birth in o['observed_births']:
        live=Path('/proc')/str(birth['pid'])/'stat'
        if live.exists(): assert live.read_text().rsplit(')',1)[1].split()[19]!=birth['saved_start_ticks']
    before=doc(p.parent/'inputs-before.json'); assert before==doc(p.parent/'inputs-after.json')
    for key,value in before.items():
        if not key.startswith('bin/'): continue
        assert pin(R/key)=={k:value[k] for k in ['bytes','sha256']},key
        if key.startswith('bin/linux-gcc/lavapipe-acqrel-ordering/'): sealed_build[key]=value
        else: products.add(key)
    owners[str(p.relative_to(R))]={'identity':pin(p),**{k:o[k] for k in ['passed','source','stop_reason','child_exit','wall_seconds','peak_sampled_rss_kib','observed_births']}}
assert owners
code=doc(W/'ordinary/ordinary-code-review.json') if (W/'ordinary/ordinary-code-review.json').exists() else None
benefit=doc(W/'readers/benefit.json') if (W/'readers/benefit.json').exists() else None
decimal=doc(W/'readers/decimal-agreement.json') if (W/'readers/decimal-agreement.json').exists() else None
all22=bool(benefit and benefit['audit_passed'] and benefit['all22_numeric_pass'])
if benefit:
    assert benefit['checks_total']==22 and decimal['agrees_on_all22'] and decimal['all22_pass']==all22
controls={str(i):({k:doc(p)[k] for k in ['passed','wall_seconds','stop_reason','child_exit','peak_sampled_rss_kib']} if (p:=W/f'runtime/control-{i}/owner.json').exists() else {'status':'not_run'}) for i in range(1,9)}
all8=(W/'readers/controls-1-8-acceptance.json').exists()
if all8: assert all22 and doc(W/'readers/controls-1-8-acceptance.json')['all8_original_controls_pass']
if not all22: assert all(controls[str(i)]=={'status':'not_run'} for i in range(1,8))
frame={'attempted':False,'accepted':False,'reason':'Original all22 and all8 gates required first.'}
p=W/'cold-frame/frame/owner.json'
if p.exists():
    assert all8
    o=doc(p); frame={k:o[k] for k in ['passed','wall_seconds','stop_reason','child_exit','peak_sampled_rss_kib']}
    frame.update(attempted=True,accepted=o['passed'],XML_present=(p.parent/'gtest.xml').exists())
findings={'candidate_source':SOURCE,'public_source_unchanged':'0109ab462edade524667d20303507aa09adcbaf6',
    'source_admission':doc(W/'review-admission.json'),'actual_provider_build':doc(W/'provider-build-result.json'),
    'actual_complete_code_review':code,'actual_new_builds':{'provider_builds':1,'Sirius_CPP_compilations':0,'shader_compilations':0,'probe_compilations':0,'frame_links':int((W/'cold-frame/build-result.json').exists())},
    'actual_mixed_host_ELF_reused':'f958a60563b64c79d3f2fb6da2815bde0effdbbe49f121034cfd212ab17b27ff',
    'original_controls':controls,'benefit_attempted':bool(benefit),'benefit_checks_passed':benefit['checks_passed'] if benefit else None,
    'all22_original_benefit_gates_pass':all22,'independent_Decimal_agrees':decimal['agrees_on_all22'] if decimal else None,
    'failed_checks':[x for x in benefit['checks']+benefit['accounting_checks'] if not x['pass']] if benefit else [],
    'whole_call_gains_percent':[{k:x[k] for k in ['direction','family','improvement_percent']} for x in benefit['checks'] if x['scope']=='complete_call_ms'] if benefit else [],
    'sole_ABBA':{'attempted':bool(benefit),'order':['a1','b1','b2','a2'],'internal_repeats_per_slot':3,'extra_retries':0},
    'all8_original_controls_admitted':all8,'cold_frame':frame,'source_or_provider_adopted_publicly':False,'full_goal_complete':False,
    'scope':'Private finite evidence only. The12s partial stack did not establish the300s timeout cause. No machine-policy selection implies numerical or performance success. Full science, July image identity, original cold frame, native-window/Radeon scenes, Mandatory and exact-main installed-release qualification remain independent obligations.'}
D.mkdir(parents=True); write(D/'findings.json',findings)
write(D/'ownership-closure.json',{'all_owned_runs_terminal_and_reaped':True,'owners_count':len(owners),'owners':owners})
diagnostics=sorted(p for p in W.rglob('*') if p.is_file())
with zipfile.ZipFile(D/'diagnostics-and-trial.zip','x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in diagnostics: z.write(p,str(p.relative_to(R)))
provider=doc(W/'provider-build-result.json'); products|={provider['library'],provider['ICD']}
for root in ['bin/linux-gcc/native-retained-ownership','bin/linux-gcc/native-cooperative-transport-intergroup']:
    products|={str(p.relative_to(R)) for p in (R/root).rglob('*') if p.is_file()}
paths={}; unique={}
for key in sorted(products):
    p=R/key; record=pin(p); record['mode']=stat.S_IMODE(p.stat().st_mode)
    if p.is_symlink(): record['symlink']=os.readlink(p)
    paths[key]=record; unique.setdefault(record['sha256'],p)
blobs={}; reused_archives={}
priors=[('08add7e5/native-original-input-ownership-finite-outcome','e4a30cfd64b6dffceb9b43690343a15e484ea8db'),
    ('bfa095ba/native-cpu-row-ownership-finite-outcome','9e01dfdca7798c79927827334798c72c2774eebc'),
    ('263abe74/lavapipe-intergroup-row-owner-finite-outcome','206ca0d4'),
    ('b1048e57/native-cooperative-transport-intergroup-finite-outcome','41101dd3')]
for directory,commit in priors:
    directory=R/'attestations/qualification'/directory; prior=doc(directory/'MANIFEST.json')
    for digest in list(unique):
        old=prior['content_addressed_blobs'].get(digest)
        if not old or old.get('reused'): continue
        archive=directory/old['archive']; key=str(archive.relative_to(R))
        if key not in reused_archives:
            assert pin(archive)==prior['archives'][old['archive']]
            b=subprocess.check_output(['git','show',commit+':'+key]); assert len(b)==archive.stat().st_size and hashlib.sha256(b).hexdigest()==pin(archive)['sha256']
            reused_archives[key]={'identity':pin(archive),'preserved_at_commit':commit}
        with zipfile.ZipFile(archive) as z:
            b=z.read(old['member']); assert len(b)==unique[digest].stat().st_size and hashlib.sha256(b).hexdigest()==digest
        blobs[digest]={'archive':key,'member':old['member'],'reused':True}; del unique[digest]
chunks=[]; chunk=[]; size=0
for digest,p in sorted(unique.items()):
    if chunk and size+p.stat().st_size>80000000: chunks.append(chunk);chunk=[];size=0
    chunk.append((digest,p)); size+=p.stat().st_size
if chunk: chunks.append(chunk)
for i,chunk in enumerate(chunks,1):
    name=f'new-actual-products-{i:02}.zip'
    with zipfile.ZipFile(D/name,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for digest,p in chunk:
            z.write(p,'sha256/'+digest); blobs[digest]={'archive':name,'member':'sha256/'+digest,'reused':False}
    assert (D/name).stat().st_size<95000000
archives={p.name:pin(p) for p in D.glob('*.zip')}
manifest={'candidate_source':SOURCE,'archives':archives,'reused_archives':reused_archives,
    'diagnostics':{str(p.relative_to(R)):pin(p) for p in diagnostics},'actual_product_paths':paths,'content_addressed_blobs':blobs,
    'provider_source_reconstruction':{'prior_bundle':'attestations/qualification/263abe74/lavapipe-intergroup-row-owner-finite-outcome/MANIFEST.json',
        'prior_manifest_identity':pin(R/'attestations/qualification/263abe74/lavapipe-intergroup-row-owner-finite-outcome/MANIFEST.json'),
        'prior_preserved_at_commit':'206ca0d4','additive_patch':pin(W/'machine-codegen.patch'),'source_contract':pin(W/'source-contract.json'),
        'instructions':'Authenticate prior bundle; reconstruct authentic Mesa25.2.8 + exact5 IEEE32 prerequisite + original intergroup prototype, then apply the additive machine-codegen patch with zero fuzz. Source contract and ROOT admission bind all8 changed files and11100-file source identity. Configured original build inputs are sealed in the build owner manifests; actual staged provider is a content-addressed product.'},
    'reconstruction':'Verify archive and blob SHA256/length, restore paths/modes then relative file symlinks. New archive names are bundle-relative; reused names are checkout-relative at their pinned pushed commits. External system/toolchain inputs are identified by original before/after manifests and are not all bundled.'}
write(D/'MANIFEST.json',manifest)
for name in archives:
    with zipfile.ZipFile(D/name) as z:
        assert z.testzip() is None
        if name=='diagnostics-and-trial.zip':
            for key,value in manifest['diagnostics'].items():
                b=z.read(key); assert len(b)==value['bytes'] and hashlib.sha256(b).hexdigest()==value['sha256']
        else:
            for key in z.namelist(): assert hashlib.sha256(z.read(key)).hexdigest()==key.split('/')[1]
assert set(blobs)=={v['sha256'] for v in paths.values()}
print(json.dumps({'bundle':str(D.relative_to(R)),'owners':len(owners),'all22':all22,'all8':all8,'cold_frame':frame['accepted'],'actual_paths':len(paths),'reused_archives':len(reused_archives),'new_archives':archives}))
