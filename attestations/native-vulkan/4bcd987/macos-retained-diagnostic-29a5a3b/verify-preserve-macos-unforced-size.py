import hashlib,json,os,pathlib,re,shutil,subprocess,tempfile,zipfile
ROOT=pathlib.Path.cwd(); scratch=ROOT/'out/macos-unforced-inline-native'
source='29a5a3b7b18256f1995b044a5cd1701fb76de9cd'
obs=scratch/'files/out/macos-retained-diagnostic/observation'
def load(p):return json.loads(p.read_text())
def ident(p):
 b=p.read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def write(p,o):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(o,indent=2)+'\n')
run=load(scratch/'run.json');job=load(scratch/'job.json');arts=load(scratch/'artifacts.json')
assert run['head_sha']==source and run['status']=='completed' and run['conclusion']=='success'
assert job['run_id']==38548451182 and job['id']==115800246512 and job['conclusion']=='success'
a=next(x for x in arts['artifacts'] if x['id']==11885297795)
assert ident(scratch/'artifact.zip')=={'bytes':8981774,'sha256':'381a836b5369edb8004f391bce5a82b9d2383e3f4740b8791f49c1b2e00ef00e'}
assert a['digest']=='sha256:'+ident(scratch/'artifact.zip')['sha256']
with zipfile.ZipFile(scratch/'artifact.zip') as z:assert z.testzip() is None
result=load(obs/'diagnostic.json');owner=load(obs/'metal-pipeline/owner.json')
assert result['runner_source']['runner_revision']==source and result['runner_source']['status']==''
assert result['status']=='metal_pipeline_observation_complete' and result['qualification_claimed'] is False and result['numerics_executed'] is False
assert result['metal_probe']['pipeline_created'] is True
assert owner['child_exit']==0 and owner['stop_reason'] is None and owner['observation_complete'] and owner['containment_complete'] and owner['retirement_census_complete']
assert not owner['remaining_observed'] and not owner['cleanup_errors'] and owner['whole_inputs_unchanged']
before=load(obs/'metal-pipeline/inputs-before.json');after=load(obs/'metal-pipeline/inputs-after.json');assert before==after
matched=0;outside=[]
prefix='/Users/runner/work/Sirius/Sirius/'
for path,expected in before.items():
 local=scratch/'files'/path[len(prefix):] if path.startswith(prefix) else None
 if local is not None and local.is_file():assert ident(local)==expected;matched+=1
 elif path.startswith(prefix+'scripts/') or path==prefix+'.github/workflows/ci.yml':
  data=subprocess.check_output(['git','show',source+':'+path[len(prefix):]])
  assert {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}==expected;matched+=1
 elif path==prefix+'out/macos-retained-diagnostic/export.zip':
  assert ident(ROOT/'attestations/native-vulkan/4bcd987/macos-retained-diagnostic-db8f439/files/out/macos-retained-diagnostic/export.zip')==expected;matched+=1
 elif path==prefix+'out/macos-retained-diagnostic/metal-source.zip':
  assert ident(ROOT/'attestations/native-vulkan/4bcd987/macos-retained-diagnostic-6e13108/artifact.zip')==expected;matched+=1
 else:outside.append({'path':path,**expected})
original=(obs/'metal-source/initialize.metal').read_bytes();derived=(obs/'metal-source/initialize-without-forced-inline.metal').read_bytes()
patt=b'static inline __attribute__((always_inline))'
assert original.count(patt)==84 and derived==original.replace(patt,b'static inline')
assert re.findall(rb'\b[0-9]+(?:\.[0-9]+)?\b',original)==re.findall(rb'\b[0-9]+(?:\.[0-9]+)?\b',derived)
events=[json.loads(x) for x in (obs/'metal-pipeline/stdout.log').read_text().splitlines()]
assert events==result['metal_probe']['events']
src=next(e for e in events if e['event']=='input');options=next(e for e in events if e['event']=='compile_options')
assert src['compiled_sha256']==ident(obs/'metal-source/initialize-without-forced-inline.metal')['sha256'] and src['removed_forced_inline_declarations']==84 and src['captured_source_modified']
assert options['optimization_level']==1 and options['language_version']==196610 and options['fast_math_enabled'] and options['preserve_invariance'] is False
assert options['changed_parameter']=={'forced_inline_declarations':{'baseline':84,'diagnostic':0}}
phases={}
for phase in ('library','function','pipeline'):
 start=next(e for e in events if e['event']=='phase_start' and e['phase']==phase)
 end=next(e for e in events if e['event']=='phase_end' and e['phase']==phase)
 assert end['success'];phases[phase]=end['monotonic_seconds']-start['monotonic_seconds']
assert events[-1]['event']=='terminal' and events[-1]['pipeline_created'] and not events[-1]['execution_claimed']
for name in ('inventory-initial','inventory-selected'):
 o=load(obs/name/'owner.json');assert o['child_exit']==0 and o['containment_complete'] and o['observation_complete'] and o['whole_inputs_unchanged']
summary={'schema':'sirius-macos-annotation-control-result-v1','source':source,'product_source':result['candidate_source'],'run':38548451182,'job':115800246512,'artifact':11885297795,'pipeline_created':True,'qualification_claimed':False,'numerics_executed':False,'native_route':'synchronous; no queue or dispatch','optimization_level':1,'removed_forced_inline_declarations':84,'phases_seconds':phases,'owner_wall_seconds':owner['wall_seconds'],'owner_peak_sampled_rss_kib':owner['peak_sampled_rss_kib'],'owner_limits':owner['limits'],'local_rechecked_bindings':matched,'remote_external_bindings_not_locally_reverified':outside,'input_manifest_equal':True,'source_bytes':ident(obs/'metal-source/initialize.metal'),'compiled_bytes':ident(obs/'metal-source/initialize-without-forced-inline.metal'),'interpretation':'Annotation-only derivative succeeds with Size on the observed Apple Paravirtual device; actual outlining, causation on other devices, MoltenVK compatibility, numerical equivalence and final release qualification are not established. A Default annotation control is a distinct next experiment.'}
write(scratch/'verification.json',summary)
bundle=ROOT/'attestations/native-vulkan/4bcd987/macos-retained-diagnostic-29a5a3b'
assert not bundle.exists();bundle.mkdir(parents=True)
for name in ('artifact.zip','artifact-metadata.json','artifacts.json','run.json','job.json','verification.json'):
 shutil.copy2(scratch/name,bundle/name)
shutil.copytree(scratch/'files',bundle/'files')
sourcepaths=['.github/workflows/ci.yml','scripts/governance/build_policy.py','scripts/diagnose-macos-retained.py','scripts/macos_metal_probe.py']
for name in sourcepaths:
 p=bundle/'source'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(subprocess.check_output(['git','show',source+':'+name]))
shutil.copy2(ROOT/'out/verify-preserve-macos-unforced-size.py',bundle/'verify-preserve-macos-unforced-size.py')
manifest={'schema':'sirius-raw-diagnostic-bundle-v1','qualification_claimed':False,'run':38548451182,'job':115800246512,'artifact':11885297795,'reused_inputs':[{'bundle':'attestations/native-vulkan/4bcd987/macos-retained-diagnostic-db8f439','path':'files/out/macos-retained-diagnostic/export.zip','sha256':'973d6ada66511b80d2dcf76e0c95a3b38367c32e08a0f16a62c20e0c0a90974d'},{'bundle':'attestations/native-vulkan/4bcd987/macos-retained-diagnostic-6e13108','path':'artifact.zip','sha256':'12060f52a9cadba5dabedb5a99b5f4293893bd066d544f2ba4e58c612d2357db'}],'files':{str(p.relative_to(bundle)):ident(p) for p in sorted(bundle.rglob('*')) if p.is_file()}}
write(bundle/'manifest.json',manifest)
# Evidence-only object tree, never checked out over product source.
index=ROOT/'out/macos-unforced-size-preservation.index';assert not index.exists()
env=dict(os.environ,GIT_INDEX_FILE=str(index));subprocess.run(['git','read-tree','--empty'],env=env,check=True)
rows=[]
for p in sorted(bundle.rglob('*')):
 if not p.is_file():continue
 sha=subprocess.check_output(['git','hash-object','-w','--',str(p)]).decode().strip()
 mode='100755' if p.stat().st_mode&0o111 else '100644'
 rows.append(f'{mode} {sha}\t{p.relative_to(ROOT)}\0')
subprocess.run(['git','update-index','-z','--index-info'],input=''.join(rows).encode(),env=env,check=True)
tree=subprocess.check_output(['git','write-tree'],env=env).decode().strip()
commit=subprocess.check_output(['git','commit-tree',tree,'-p',source],input=b'Preserve successful annotation-only native Metal Size control (#48)\n',env=env).decode().strip()
tag='evidence/issue-48-macos-diagnostic-29a5a3b';subprocess.run(['git','tag',tag,commit],check=True)
index.unlink()
print(json.dumps({'summary':summary,'bundle':str(bundle.relative_to(ROOT)),'files':len(manifest['files'])+1,'manifest':ident(bundle/'manifest.json'),'tag':tag,'commit':commit},indent=2))
