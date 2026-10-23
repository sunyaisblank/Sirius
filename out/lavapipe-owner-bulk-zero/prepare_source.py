"""Restore optimized intergroup basis and express only its large zero store as memset."""
from pathlib import Path
import difflib, hashlib, json, subprocess, zipfile
R=Path.cwd();W=Path(__file__).resolve().parent
S=R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
OLD=R/'out/lavapipe-owner-machine-codegen'
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
def doc(p):return json.loads(p.read_text())
review=doc(OLD/'review-admission.json'); assert review['passed']
current={str(p.relative_to(S)):pin(p) for p in S.rglob('*') if p.is_file()}
b=json.dumps(current,sort_keys=True,separators=(',',':')).encode()
assert {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}==review['full_current_source_identity']
assert pin(OLD/'machine-codegen.patch')==review['machine_codegen_patch']
assert not (W/'bulk-zero.patch').exists()
subprocess.run(['patch','--reverse','--batch','--fuzz=0','-p1','-d',str(S),'-i',str(OLD/'machine-codegen.patch')],check=True,capture_output=True)
archive=R/review['basis']['archive'];assert pin(archive)=={k:review['basis'][k] for k in ['bytes','sha256']}
with zipfile.ZipFile(archive) as z:before=json.loads(z.read(review['basis']['member']))
prefix=str(S.relative_to(R))+'/'
prior={k[len(prefix):]:{j:v[j] for j in ['bytes','sha256']} for k,v in before.items() if k.startswith(prefix)}
assert len(prior)==11100
assert prior=={str(p.relative_to(S)):pin(p) for p in S.rglob('*') if p.is_file()}
files=['src/gallium/auxiliary/gallivm/lp_bld_nir_soa.c','src/gallium/frontends/lavapipe/lvp_pipeline.c','src/gallium/frontends/lavapipe/lvp_device.c','src/gallium/drivers/llvmpipe/lp_screen.c']
basis={k:(S/k).read_text() for k in files}
p=S/files[0];source=p.read_text()
old='''      LLVMValueRef reg_alloc = lp_build_alloca(bld.base.gallivm,
                                               type, "reg");'''
new='''      LLVMValueRef reg_alloc;
      bool bulk_zero = false;
#if LP_VULKAN_IEEE32_ENABLED && !GALLIVM_USE_ORCJIT && LLVM_VERSION_MAJOR == 21
      /* Exact admitted row owners retain all private scratch words, including
       * padding. Express the same all-bits-zero store as one bulk operation;
       * do not make the optimizer and machine allocator process thousands of
       * independent vector-zero stores. Ordinary variants keep their lowering.
       */
      if (shader->info.stage == MESA_SHADER_COMPUTE &&
          shader->info.cs.llvmpipe_intergroup_row_owner &&
          (shader->options->divergence_analysis_options &
           nir_divergence_multiple_workgroup_per_compute_subgroup) &&
          LLVMGetTypeKind(type) == LLVMArrayTypeKind) {
         LLVMTypeRef element = LLVMGetElementType(type);
         unsigned count = LLVMGetArrayLength(type);
         bulk_zero = (count == 3060 || count == 3145) &&
                     LLVMGetTypeKind(element) == LLVMVectorTypeKind &&
                     LLVMGetVectorSize(element) == 8 &&
                     LLVMGetTypeKind(LLVMGetElementType(element)) == LLVMIntegerTypeKind &&
                     LLVMGetIntTypeWidth(LLVMGetElementType(element)) == 32;
      }
#endif
      if (bulk_zero) {
         /* Same entry-block allocation, original insertion point and complete
          * ABI allocation size. No rounding, skipped byte or lifetime change.
          */
         reg_alloc = lp_build_alloca_undef(gallivm, type, "reg");
         LLVMBuildMemSet(gallivm->builder, reg_alloc,
                         LLVMConstInt(LLVMInt8TypeInContext(gallivm->context), 0, false),
                         LLVMConstInt(LLVMInt64TypeInContext(gallivm->context),
                                      LLVMABISizeOfType(gallivm->target, type), false),
                         LLVMABIAlignmentOfType(gallivm->target, type));
      } else {
         reg_alloc = lp_build_alloca(gallivm, type, "reg");
      }'''
assert source.count(old)==1; p.write_text(source.replace(old,new))
p=S/files[1];source=p.read_text()
line='      {0xc0,0xf3,0x8b,0x8b,0x39,0xd1,0x21,0x36,0x21,0x6a,0xf9,0x8f,0x34,0x65,0xf3,0x46,0x46,0x30,0x06,0x28},\n'
assert source.count(line)==1;p.write_text(source.replace(line,'').replace('These three complete-row shaders','These two complete-row shaders'))
p=S/files[2];source=p.read_text();assert source.count('SIR-IWG-25.2.8-1')==1;p.write_text(source.replace('SIR-IWG-25.2.8-1','SIR-IWG-25.2.8-3'))
p=S/files[3];source=p.read_text();assert source.count('Sirius-intergroup-owner-v1')==1;p.write_text(source.replace('Sirius-intergroup-owner-v1','Sirius-intergroup-owner-bulkzero-v3'))
patch=''.join(''.join(difflib.unified_diff(basis[k].splitlines(True),(S/k).read_text().splitlines(True),fromfile='a/'+k,tofile='b/'+k)) for k in files)
(W/'bulk-zero.patch').write_text(patch)
contract={'scope':'Exact E/DP intergroup-owner private uint array initialization only; restored original Default optimized machine pipeline. All zero bytes retained.',
 'source_basis':{'prior_intergroup_ROOT_review':{'path':'out/lavapipe-intergroup-row-owner/review-admission.json',**pin(R/'out/lavapipe-intergroup-row-owner/review-admission.json')},
 'prior_intergroup_contract':{'path':'out/lavapipe-intergroup-row-owner/source-contract.json',**pin(R/'out/lavapipe-intergroup-row-owner/source-contract.json')},
 'prior_intergroup_patch':{'path':'out/lavapipe-intergroup-row-owner/prototype.patch',**pin(R/'out/lavapipe-intergroup-row-owner/prototype.patch')},
 'IEEE32_prerequisite':{'path':'out/lavapipe-acqrel-ordering/ieee32-prerequisite.patch',**pin(R/'out/lavapipe-acqrel-ordering/ieee32-prerequisite.patch')}},
 'removed_rejected_None_policy':{'patch':{'path':str((OLD/'machine-codegen.patch').relative_to(R)),**pin(OLD/'machine-codegen.patch')},'complete_11100_file_prior_basis_restored_before_new_edit':True},
 'modified_files':{k:{'basis':prior[k],'candidate':{'path':str((S/k).relative_to(R)),**pin(S/k)}} for k in files},
 'prototype_patch':pin(W/'bulk-zero.patch'),
 'actual_intended_host':doc(OLD/'source-contract.json')['actual_intended_host'],
 'zero_extent':{'Endpoint':{'vectors':3060,'vector_bytes':32,'bytes':97920},'DopriPhase':{'vectors':3145,'vector_bytes':32,'bytes':100640}},
 'precondition':'Exact E/DP SPIR-V and existing structural/dispatch gates, LLVM21/MCJIT0/IEEE32, multiple-workgroup NIR option, one flat array of3060or3145 vectors of8xi32.',
 'postcondition':'Identical complete bitwise-zero array before the original first use; same entry allocation, builder insertion point, extent/alignment and original dataflow. No omitted initialization or undefined payload assumption.',
 'frame_conditions':['No shader/host/module/ABI/group/field/arithmetic/mode/count/admission/refusal/rollback change.','No shared/global optimization mutation; original LLVM pass pipeline and CodeGenDefault restored.','Ordinary key0 and both Transport hashes use original allocation+aggregate-null store.','Private UUID3 and disk namespace bulkzero-v3 isolate machine objects.'],
 'actual_code_gate':'Before C8 require selected actual optimized IR one complete memset for each E/DP scratch array with no expanded prologue-zero stores; compare remaining IR by a narrowly declared equivalent-clear normalization. Actual called Default machine bytes/clear lowering and ordinary T isolation/features/FP-state must be bound. Runtime benefit not predicted.',
 'builds_or_workloads':0,'qualified':False}
(W/'source-contract.json').write_text(json.dumps(contract,indent=2)+'\n')
print(json.dumps({'restored_prior_files':11100,'changed_files':files,'patch':pin(W/'bulk-zero.patch'),'shader_or_host_changes':False,'builds_or_workloads':0}))
