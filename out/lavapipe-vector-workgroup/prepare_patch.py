from pathlib import Path
import difflib,hashlib,json,re
R=Path.cwd();W=Path(__file__).resolve().parent;S=R/'bin/toolchains/mesa-25.2.8'
NAMES=['src/gallium/drivers/llvmpipe/lp_state_cs.c','src/gallium/auxiliary/gallivm/lp_bld_nir_soa.c','src/gallium/auxiliary/gallivm/lp_bld_tgsi.h']
original={n:(S/n).read_text() for n in NAMES};modified=original.copy()
helper=r'''
/* A complete static workgroup may share one LLVM vector even when the
 * advertised subgroup width stays eight. No subgroup-sensitive shader may
 * enter this path. General, variable, mesh and unsupported CPU workgroups
 * retain the original coroutine implementation.
 * This private candidate is disabled in the matched baseline build. */
#ifndef LP_SINGLE_VECTOR_WORKGROUP
#define LP_SINGLE_VECTOR_WORKGROUP 0
#endif
static bool
single_vector_workgroup(const nir_shader *nir)
{
   if (!LP_SINGLE_VECTOR_WORKGROUP ||
       nir->info.stage != MESA_SHADER_COMPUTE ||
       !nir->info.uses_control_barrier || nir->info.workgroup_size_variable ||
       nir->info.workgroup_size[0] != 16 ||
       nir->info.workgroup_size[1] != 1 || nir->info.workgroup_size[2] != 1 ||
       nir->info.shared_size > 16384 || lp_native_vector_width != 256 ||
       util_get_cpu_caps()->max_vector_bits < 512 ||
       (nir->info.bit_sizes_float & ~32u) ||
       (nir->info.bit_sizes_int & ~(1u | 32u)) ||
       nir->info.uses_wide_subgroup_intrinsics)
      return false;

   nir_foreach_function_impl(impl, nir) {
      nir_foreach_block(block, impl) {
         nir_foreach_instr(instr, block) {
            if (instr->type != nir_instr_type_intrinsic)
               continue;
            const nir_intrinsic_instr *intr = nir_instr_as_intrinsic(instr);
            if (nir_intrinsic_has_semantic(intr, NIR_INTRINSIC_SUBGROUP) ||
                strstr(nir_intrinsic_infos[intr->intrinsic].name, "subgroup"))
               return false;
            if (intr->intrinsic == nir_intrinsic_barrier &&
                nir_intrinsic_execution_scope(intr) != SCOPE_NONE &&
                nir_intrinsic_execution_scope(intr) != SCOPE_WORKGROUP)
               return false;
         }
      }
   }
   return true;
}

'''
n=NAMES[0];s=modified[n]
s=s.replace('#include "util/u_memory.h"','#include "util/u_memory.h"\n#include "util/u_cpu_detect.h"',1)
mark='static void\ngenerate_compute('
assert s.count(mark)==1;s=s.replace(mark,helper+mark,1)
old='   bool use_coro = nir->info.uses_control_barrier || is_mesh;';assert old in s
s=s.replace(old,'   const bool single_vector = single_vector_workgroup(nir);\n   bool use_coro = (nir->info.uses_control_barrier && !single_vector) || is_mesh;',1)
old='   cs_type.length = MIN2(lp_native_vector_width / 32, 16); /* n*4 elements per vector */';assert old in s
s=s.replace(old,'   cs_type.length = single_vector ? 16 : MIN2(lp_native_vector_width / 32, 16);',1)
# Both the entry point and any called function get the same explicit barrier
# lowering contract. Other seven original parameter producers start zeroed.
assert s.count('params.type = cs_type;')==2
s=s.replace('params.type = cs_type;','params.type = cs_type;\n         params.single_vector_workgroup = single_vector;',1)
i=s.index('params.type = cs_type;',s.index('params.single_vector_workgroup'))
s=s[:i]+s[i:].replace('params.type = cs_type;','params.type = cs_type;\n      params.single_vector_workgroup = single_vector;',1)
modified[n]=s
n=NAMES[1];s=modified[n]
old='   const struct lp_build_coro_suspend_info *coro;';assert s.count(old)==1
s=s.replace(old,old+'\n   bool single_vector_workgroup;',1)
old='   if (exec_scope != SCOPE_NONE) {';assert s.count(old)==1
s=s.replace(old,'''   /* One LLVM vector executes every workgroup invocation before advancing
    * beyond this point. Keep the original memory fence above; only the
    * between-vector suspension is unnecessary. The eligibility predicate
    * excludes subgroup scope and every subgroup-observable shader operation. */
   assert(!bld->single_vector_workgroup ||
          exec_scope == SCOPE_NONE || exec_scope == SCOPE_WORKGROUP);
   if (exec_scope != SCOPE_NONE && !bld->single_vector_workgroup) {''',1)
old='   bld.coro = params->coro;';assert old in s;s=s.replace(old,old+'\n   bld.single_vector_workgroup = params->single_vector_workgroup;',1)
modified[n]=s
n=NAMES[2];s=modified[n]
start=s.index('struct lp_build_tgsi_params {');end=s.index('\n};',start)
s=s[:end]+'\n   /* Compute-only: the complete static workgroup occupies one LLVM vector. */\n   bool single_vector_workgroup;'+s[end:]
modified[n]=s
records={};patch=[]
for n in NAMES:
 for directory,value in [('original',original[n]),('candidate',modified[n])]:
  p=W/directory/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(value)
 patch.extend(difflib.unified_diff(original[n].splitlines(keepends=True),modified[n].splitlines(keepends=True),fromfile='a/'+n,tofile='b/'+n))
 records[n]={arm:{'bytes':len(value.encode()),'sha256':hashlib.sha256(value.encode()).hexdigest()} for arm,value in [('original',original[n]),('candidate',modified[n])]}
(W/'provider.patch').write_text(''.join(patch))
zeroed=[]
for p in (S/'src').rglob('*.c'):
 text=p.read_text(errors='replace')
 for m in re.finditer(r'struct lp_build_tgsi_params\s+(\w+)\s*;',text):
  fragment=text[m.end():m.end()+200]
  assert re.match(r'\s*memset\(&'+m[1]+r', 0, sizeof\('+m[1]+r'\)\);',fragment)
  zeroed.append({'file':str(p.relative_to(S)),'line':text[:m.start()].count('\n')+1})
assert len(zeroed)==7
(W/'source-contract-review.json').write_text(json.dumps({'schema':'sirius-lavapipe-single-vector-source-contract-v1','source_screen_passed':True,'not_a_machine_checked_proof':True,'Sirius_source_changed':False,'shader_or_program_or_arithmetic_changed':False,'provider_patch_adopted':False,'changed_vendor_sources':records,'parameter_producers_zero_initialised':zeroed,'eligibility':'Compute; constant16x1x1; original control barrier; shared<=16KiB; defaultnative256; CPUmax>=512; no subgroup semantic flag/builtin; only NONE/WORKGROUP execution barriers. Unsupported cases preserve original path.','execution_contract':'All16 invocations run in one LLVM vector; original structured control/masking reconverges the workgroup before the original memory fence; execution suspension between vectors is unnecessary. Every memory fence, memory scope/semantics, atomic, shared/global access, NIR barrier and program remains unchanged. Verify ordinary called LLVM and original complete readbacks before runtime/cost acceptance.','NIR_ownership':'No NIR clone or mutation; explicit zero-default lowering context flag only, set in both compute entry/called-function producers.','remaining_uncertainty':'Actual supported build, executed vector/barrier lowering, complete numerical/refusal/rollback and original22 benefit remain unproved.'},indent=2)+'\n')
print(json.dumps({'changed_vendor_files':len(NAMES),'new_lines':sum(len(modified[n].splitlines())-len(original[n].splitlines()) for n in NAMES),'zeroed_parameter_producers':len(zeroed),'no_Sirius_runtime_or_shader_changes':True}))
