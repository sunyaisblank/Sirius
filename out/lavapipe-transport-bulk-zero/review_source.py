"""Authenticate the exact source delta and full unchanged source complement."""
from pathlib import Path
import hashlib,json,shutil,subprocess,tempfile
R=Path.cwd();W=Path(__file__).resolve().parent;S=R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
def doc(p):return json.loads(p.read_text())
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
def canonical(d):
 b=json.dumps(d,sort_keys=True,separators=(',',':')).encode()
 return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
c=doc(W/'source-contract.json');a=c['basis'];p=R/a['review_path']
assert pin(p)==a['review_identity'] and doc(p)['full_current_source_identity']==a['full_source_identity']
assert c['prototype_patch']==pin(W/'bulk-zero.patch')
current={str(p.relative_to(S)):pin(p) for p in S.rglob('*') if p.is_file()};assert len(current)==11100
basis=dict(current)
for k,item in c['modified_files'].items():
 assert current[k]=={j:item['candidate'][j] for j in ['bytes','sha256']};basis[k]=item['basis']
assert canonical(basis)==a['full_source_identity']
with tempfile.TemporaryDirectory(dir=W,prefix='reverse-check-') as folder:
 folder=Path(folder)
 for k in c['modified_files']:
  p=folder/k;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(S/k,p)
 subprocess.run(['patch','--reverse','--batch','--fuzz=0','-p1','-d',str(folder),'-i',str(W/'bulk-zero.patch')],check=True,capture_output=True)
 assert not list(folder.rglob('*.rej'))
 for k in c['modified_files']:assert pin(folder/k)==basis[k]
soa=(S/'src/gallium/auxiliary/gallivm/lp_bld_nir_soa.c').read_text()
assert soa.count('(count == 2295 || count == 3060 || count == 3145)')==1
for s in ['LLVMBuildMemSet','LLVMABISizeOfType(gallivm->target, type)','LLVMABIAlignmentOfType(gallivm->target, type)','LP_VULKAN_IEEE32_ENABLED && !GALLIVM_USE_ORCJIT && LLVM_VERSION_MAJOR == 21']:
 assert s in soa
t=c['Transport'];assert t['vectors']*t['lanes']*t['bits']//8==t['bytes']==73440
for k in ['src/gallium/auxiliary/gallivm/lp_bld_init.h','src/gallium/auxiliary/gallivm/lp_bld_init.c','src/gallium/drivers/llvmpipe/lp_state_cs.c','src/gallium/auxiliary/gallivm/lp_bld_misc.cpp']:
 assert 'codegen_none_after_ir' not in (S/k).read_text()
result={'passed':True,'prototype_patch':pin(W/'bulk-zero.patch'),'source_contract':pin(W/'source-contract.json'),
 'source_files':{k:current[k] for k in sorted(c['modified_files'])},'full_current_source_identity':canonical(current),
 'basis':a,'unchanged_remaining_files':11096,'all11100_source_files_authenticated':True,
 'independent_source_review':'Existing divide_contract_review independently grounded exact fullzero Transport candidate; ROOT inspected literal4-file delta. No runtime conclusion.',
 'reviewed_contract':c['postcondition'],'required_code_gate':c['required_code_gate'],'runtime_or_benefit_admitted':False}
p=W/'review-admission.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'source_review_passed':True,'sealed_source_files':11100,'changed_files':4,'full_zero_bytes':73440,'runtime_admitted':False}))
