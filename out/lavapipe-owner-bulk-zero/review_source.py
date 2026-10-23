"""Authenticate all source and the equivalent complete-clear contract before build."""
from pathlib import Path
import hashlib,json,subprocess,tempfile,shutil,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent;S=R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
def doc(p):return json.loads(p.read_text())
def data_pin(b):return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def pin(p):return data_pin(p.read_bytes())
def identity(d):return {k:d[k] for k in ['bytes','sha256']}
c=doc(W/'source-contract.json');assert pin(W/'bulk-zero.patch')==c['prototype_patch']
for item in c['source_basis'].values():assert pin(R/item['path'])==identity(item)
old=doc(R/'out/lavapipe-owner-machine-codegen/review-admission.json');a=old['basis'];archive=R/a['archive']
assert pin(archive)==identity(a)
with zipfile.ZipFile(archive) as z:b=z.read(a['member'])
assert data_pin(b)==a['member_identity']
prefix=str(S.relative_to(R))+'/'
prior={k[len(prefix):]:identity(v) for k,v in json.loads(b).items() if k.startswith(prefix)}
current={str(p.relative_to(S)):pin(p) for p in S.rglob('*') if p.is_file()}
assert len(prior)==len(current)==11100 and prior.keys()==current.keys()
changed={k for k in current if current[k]!=prior[k]};assert changed==set(c['modified_files']) and len(changed)==4
for k,item in c['modified_files'].items():assert prior[k]==identity(item['basis']) and current[k]==identity(item['candidate']),k
with tempfile.TemporaryDirectory(dir=W,prefix='reverse-check-') as temporary:
 temporary=Path(temporary)
 for k in changed:
  p=temporary/k;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(S/k,p)
 subprocess.run(['patch','--reverse','--batch','--fuzz=0','-p1','-d',str(temporary),'-i',str(W/'bulk-zero.patch')],check=True,capture_output=True)
 assert not list(temporary.rglob('*.rej'))
 for k in changed:assert pin(temporary/k)==prior[k],k
for k in old['source_files']:
 if k not in changed:assert current[k]==prior[k],k
assert all('intergroup_codegen_none_after_ir' not in (S/k).read_text() for k in ['src/gallium/auxiliary/gallivm/lp_bld_init.h','src/gallium/auxiliary/gallivm/lp_bld_init.c','src/gallium/drivers/llvmpipe/lp_state_cs.c'])
assert 'lp_set_mcjit_codegen_none_after_ir' not in (S/'src/gallium/auxiliary/gallivm/lp_bld_misc.cpp').read_text()
for shape in c['zero_extent'].values():assert shape['vectors']*8*4==shape['bytes'] and shape['vector_bytes']==32
commands=doc(R/'bin/linux-gcc/lavapipe-acqrel-ordering/compile_commands.json')
rows=[r for r in commands if Path(r['file']).name=='lp_bld_nir_soa.c'];assert len(rows)==1
assert all(flag in rows[0]['command'] for flag in ['-DGALLIVM_USE_ORCJIT=0','-DLP_VULKAN_IEEE32=1'])
host=c['actual_intended_host'];record=host['matched_result'];assert pin(R/record['path'])==identity(record)
matched=doc(R/record['path']);exe=R/matched['executable'];assert pin(exe)['sha256']==host['executable_sha256']==matched['executable_sha256']
for kind,item in host['modules'].items():
 module=R/item['module'];assert pin(module)==identity(item) and exe.read_bytes().count(module.read_bytes())==1,kind
result={'passed':True,'prototype_patch':pin(W/'bulk-zero.patch'),'source_contract':pin(W/'source-contract.json'),
 'source_files':{k:current[k] for k in sorted(changed)},'full_current_source_identity':data_pin(json.dumps(current,sort_keys=True,separators=(',',':')).encode()),
 'basis':a,'changed_files':sorted(changed),'unchanged_remaining_files':11100-4,
 'rejected_None_policy_completely_removed':True,'original_Default_machine_pipeline_and_IR_pass_calls_restored':True,
 'reviewed_contract':'Exact two admitted compute variants and exact array shape only. Same entry allocation, original current-builder store point, all ABI bytes zero/nonvolatile, target-derived alignment. No lifetime/mask/field/math/ABI/fallback change. All other callers retain original aggregate-null initialization.',
 'independent_source_review':'divide_contract_review reviewed source Hoare contract and concrete patch before build; actual optimizer retention/publication/runtime gates still required.',
 'actual_code_gate':c['actual_code_gate'],'new_Sirius_CPP_or_shader_compilations':0,'runtime_or_benefit_admitted':False}
p=W/'review-admission.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'source_review_passed':True,'sealed_source_files':11100,'changed_files':4,'full_zero_bytes':c['zero_extent'],'runtime_admitted':False}))
