"""Extend the authenticated bulk-zero source to the exact full Transport row owner."""
from pathlib import Path
import difflib,hashlib,json
R=Path.cwd();W=Path(__file__).resolve().parent
S=R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
OLD=R/'out/lavapipe-owner-bulk-zero'
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
def canonical(d):
 b=json.dumps(d,sort_keys=True,separators=(',',':')).encode()
 return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
review=json.loads((OLD/'review-admission.json').read_text());assert review['passed']
current={str(p.relative_to(S)):pin(p) for p in S.rglob('*') if p.is_file()}
assert len(current)==11100 and canonical(current)==review['full_current_source_identity']
assert not (W/'bulk-zero.patch').exists()
files=['src/gallium/auxiliary/gallivm/lp_bld_nir_soa.c','src/gallium/frontends/lavapipe/lvp_pipeline.c','src/gallium/frontends/lavapipe/lvp_device.c','src/gallium/drivers/llvmpipe/lp_screen.c']
basis={k:(S/k).read_text() for k in files}
p=S/files[0];s=p.read_text();old='(count == 3060 || count == 3145)'
assert s.count(old)==1;p.write_text(s.replace(old,'(count == 2295 || count == 3060 || count == 3145)'))
p=S/files[1];s=p.read_text();old='      {0xe4,0xf3,0x2b,0xfe,0x7d,0xb7,0x3d,0xad,0x21,0xd7,0x35,0x86,0xf1,0xbf,0x77,0x17,0x2a,0xa3,0x90,0xca},'
assert s.count(old)==1
new='      {0xc0,0xf3,0x8b,0x8b,0x39,0xd1,0x21,0x36,0x21,0x6a,0xf9,0x8f,0x34,0x65,0xf3,0x46,0x46,0x30,0x06,0x28},\n'+old
p.write_text(s.replace(old,new).replace('These two complete-row shaders','These three complete-row shaders'))
for k,old,new in [(files[2],'SIR-IWG-25.2.8-3','SIR-IWG-25.2.8-4'),(files[3],'Sirius-intergroup-owner-bulkzero-v3','Sirius-intergroup-transport-bulkzero-v4')]:
 p=S/k;s=p.read_text();assert s.count(old)==1;p.write_text(s.replace(old,new))
patch=''.join(''.join(difflib.unified_diff(basis[k].splitlines(True),(S/k).read_text().splitlines(True),fromfile='a/'+k,tofile='b/'+k)) for k in files)
(W/'bulk-zero.patch').write_text(patch)
c={'scope':'Exact complete row-Transport initialization extension; previous E/DP implementation unchanged. No predicted numerical, benefit or frame result.',
 'basis':{'review_path':str((OLD/'review-admission.json').relative_to(R)),'review_identity':pin(OLD/'review-admission.json'),
 'full_source_identity':review['full_current_source_identity'],
 'preserved_at_commit':'48adb23d','manifest':'attestations/qualification/337d0bae/lavapipe-owner-bulk-zero-finite-outcome/MANIFEST.json'},
 'prototype_patch':pin(W/'bulk-zero.patch'),
 'modified_files':{k:{'basis':current[k],'candidate':{'path':str((S/k).relative_to(R)),**pin(S/k)}} for k in files},
 'Transport':{'SPIRV_SHA1':'c0f38b8b39d12136216af98f3465f34646300628','vectors':2295,'lanes':8,'bits':32,'bytes':73440,'ABI_alignment':32},
 'precondition':'Exact row-T hash plus original one-function/barrierless structural gate, compute/marker/multiple-workgroup divergence/IEEE32/LLVM21/MCJIT0; flat2295x8xi32 shape.',
 'postcondition':'Same entry alloca, original builder initialization point and complete ABI extent/alignment; every active/inactive/padding byte zero before original use; original full arithmetic/masks/IDs/lifetime remain.',
 'frame_conditions':['Original Default IR and machine pipeline, all40fields and all46inputs, consumers/precision/refusal/rollback/caches/governors/counts unchanged.','No shader or Sirius C++ changes or compilations; archived complete all-row host can be authenticated and restored if code gate passes.','Previous rejected18/22 Transport source did not contain this bulk-zero representation; no unchanged retry.','Fresh private UUID4 and transport-bulkzero-v4 namespace prevent stale machine objects.'],
 'required_code_gate':'One actual4-byte guard capture must bind complete73440-byte clear and allzero return bytes, resolved libc publication, same remaining optimized IR/attributes except explicitly justified defined Boolean-mask identities, actual Default engine+greedy, original module/features/properties/FP restore. No numerical gate until this passes.',
 'independent_premise':'Saved exact T-row optimized IR expands2295zero-vector stores; actual oldjit has no memset relocation and expanded store prologue. This justifies a code screen, without predicting runtime benefit.',
 'new_Sirius_CPP_or_shader_compilations':0,'builds_or_workloads':0,'qualified':False}
(W/'source-contract.json').write_text(json.dumps(c,indent=2)+'\n')
print(json.dumps({'authenticated_basis_files':len(current),'changed_files':files,'full_zero_bytes':73440,'builds_or_workloads':0}))
