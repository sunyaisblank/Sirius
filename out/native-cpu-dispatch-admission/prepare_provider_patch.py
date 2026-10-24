"""Emit complete-dispatch admission before LLVM, with original refused-job path.

No shared provider source mutation, build, provider loading or workload occurs.
"""
from pathlib import Path
import difflib
import hashlib
import json

root = Path.cwd()
work = Path(__file__).resolve().parent
provider = work / 'provider'
basis = root / 'out/native-cpu-workgroup-scheduler'
source_root = root / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
provider.mkdir(exist_ok=True)


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def canonical(records):
    return pin(json.dumps(records, sort_keys=True, separators=(',', ':')).encode())


def replace(text, before, after):
    assert text.count(before) == 1, before
    return text.replace(before, after)


def function(text, marker):
    start = text.index(marker)
    opening = text.index('{', start)
    depth, end = 1, opening + 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}')
        end += 1
    return start, end


assert not (provider / 'provider.patch').exists()
original = json.loads((basis / 'provider/source-contract.json').read_text())
before = {str(p.relative_to(source_root)): pin(p.read_bytes())
          for p in source_root.rglob('*') if p.is_file()}
assert canonical(before) == original['full_candidate_source_identity']
compiled = json.loads((work / 'batch-compile-result.json').read_text())
owner = json.loads((work / 'build/compile/owner.json').read_text())
assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert compiled['all_prepared_source_bytes_unchanged']
assert len(compiled['objects']) == 1 and compiled['objects'][0]['exit'] == 0
objects = original['native_objects'] + [compiled['objects'][0]['object']]
assert len(objects) == 4
for item in objects:
    assert pin((root / item['path']).read_bytes()) == {k: item[k] for k in ('bytes', 'sha256')}

cs_name = 'src/gallium/drivers/llvmpipe/lp_state_cs.c'
cs_old = (source_root / cs_name).read_text()
cs = replace(cs_old, '   bool intergroup_row_owner;\n',
             '   bool intergroup_row_owner;\n   bool native_preflight_passed;\n')
start, end = function(cs, 'static void\nllvmpipe_cs_update_derived(')
derived = cs[start:end]
variant = derived.index('   if (llvmpipe->cs_dirty & (LP_CSNEW_CS |')
resources = derived[:variant].replace('llvmpipe_cs_update_derived(', 'llvmpipe_cs_update_resources(', 1)
resources += '''   /* Resource flags are consumed; preserve any pending variant update. */
   llvmpipe->cs_dirty = (llvmpipe->cs_dirty & (LP_CSNEW_CS |
                             LP_CSNEW_IMAGES |
                             LP_CSNEW_SAMPLER_VIEW |
                             LP_CSNEW_SAMPLER)) ? LP_CSNEW_CS : 0;
}


static void
llvmpipe_cs_update_derived(struct llvmpipe_context *llvmpipe)
{
   llvmpipe_cs_update_resources(llvmpipe);
   if (llvmpipe->cs_dirty & LP_CSNEW_CS)
      llvmpipe_update_cs(llvmpipe);
   llvmpipe->cs_dirty = 0;
}'''
cs = cs[:start] + resources + cs[end:]
cs = replace(cs, '''   if (!native_handled)
      variant->jit_function''', '''   /* Complete immutable admission proves every scheduled packet is handled.
    * Arithmetic refusal is handled execution. A false here would violate
    * that internal lifetime/guard invariant; never dereference a NULL JIT.
    */
   if (job_info->native_preflight_passed && !native_handled)
      abort();
   if (!native_handled)
      variant->jit_function''')
cs = replace(cs, '''   /* Native scalar groups need one original workgroup per pool task.
    * Keep the existing ordinary-ID variant for every packet refusal.
    */
''', '')
cs = replace(cs, '      !llvmpipe->cs->base.ir.nir->info.cs.llvmpipe_native_cpu_stage &&\n', '')
native_start = cs.index('#if LP_VULKAN_IEEE32_ENABLED\n   const nir_shader *native_nir')
native_end = cs.index('#endif', native_start) + len('#endif')
native_block = cs[native_start:native_end]
coarse = replace(native_block, '   if (native_nir->info.stage', '   native_direct = (native_nir->info.stage')
coarse = replace(coarse, ')\n      job_info.native_cpu_stage = native_nir->info.cs.llvmpipe_native_cpu_stage;', ');')
coarse = coarse.replace('job_info.grid_size', 'info->grid')
cs = cs[:native_start] + '''#if LP_VULKAN_IEEE32_ENABLED
   if (native_direct) {
      job_info.native_cpu_stage = native_nir->info.cs.llvmpipe_native_cpu_stage;
      if (lp_native_cpu_preflight(&job_info)) {
         job_info.native_preflight_passed = true;
         job_info.current->variant = NULL;
         /* No variant is bound. A later refusal must select one even when
          * resource bindings and the original fallback key are unchanged.
          */
         llvmpipe->cs_dirty |= LP_CSNEW_CS;
         intergroup = false;
      } else {
         /* The original eligible key and JIT exist before any output write. */
         llvmpipe_cs_update_derived(llvmpipe);
      }
   }
#endif''' + cs[native_end:]
cs = replace(cs, '   llvmpipe_cs_update_derived(llvmpipe);\n\n   fill_grid_size',
             '   bool native_direct = false;\n' + coarse + '''
   /* Preserve derived/JIT/grid resolution order for all other dispatches. */
   if (native_direct)
      llvmpipe_cs_update_resources(llvmpipe);
   else
      llvmpipe_cs_update_derived(llvmpipe);

   fill_grid_size''')
assert cs.count('llvmpipe_cs_update_derived(llvmpipe);') == 2
assert '+= (uint64_t)num_tasks * info->block[0]' in cs

bridge_name = 'src/gallium/drivers/llvmpipe/lp_native_cpu_bridge.h'
bridge_old = (source_root / bridge_name).read_text()
bridge = replace(bridge_old,
                 'lp_native_cpu_packet(const struct lp_cs_job_info *job,\n                     uint32_t x, uint32_t y, uint32_t z)',
                 'lp_native_cpu_packet_mode(const struct lp_cs_job_info *job,\n                          uint32_t x, uint32_t y, uint32_t z, bool preflight)')
switch = bridge[bridge.index('   switch (job->native_cpu_stage) {'):bridge.rindex('   return false;')]
pure_switch = switch.replace('_try_packet(', '_admit_packet(')
bridge = replace(bridge, switch, '   if (preflight) {\n' + pure_switch + '   }\n' + switch)
bridge += '''
static bool
lp_native_cpu_packet(const struct lp_cs_job_info *job,
                     uint32_t x, uint32_t y, uint32_t z)
{
   return lp_native_cpu_packet_mode(job, x, y, z, false);
}

/* All job/descriptor/header/table/extent/alias checks are shared with the
 * worker. At x=0/groups=1, PASS proves every original 1D workgroup packet
 * fits the same immutable capacity. This function writes no output.
 */
static bool
lp_native_cpu_preflight(const struct lp_cs_job_info *job)
{
   return lp_native_cpu_packet_mode(job, 0, 0, 0, true);
}
'''
abi_name = 'src/gallium/drivers/llvmpipe/sirius_native_bridge.h'
abi_old = (source_root / abi_name).read_text()
exports = '\n// Pure full-dispatch admission: no floating arithmetic or output writes.\n'
for slug in ('transport', 'endpoint', 'dopri_phase'):
    exports += f'bool sirius_native_{slug}_admit_packet(SIRIUS_NATIVE_PACKET_ARGUMENTS) SIRIUS_NATIVE_NOEXCEPT;\n'
abi = replace(abi_old, '\n#undef SIRIUS_NATIVE_PACKET_ARGUMENTS', exports + '\n#undef SIRIUS_NATIVE_PACKET_ARGUMENTS')
meson_name = 'src/gallium/targets/lavapipe/meson.build'
meson_old = (source_root / meson_name).read_text()
old_last = "    '" + str(root / objects[2]['path']) + "'\n"
new_last = old_last.rstrip('\n') + ",\n    '" + str(root / objects[3]['path']) + "'\n"
meson = replace(meson_old, old_last, new_last)
new_text = {cs_name: cs, bridge_name: bridge, abi_name: abi, meson_name: meson}
for name, a, b in [
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8E', 'SIR-CPU-25.2.8F'),
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-workgroup-scheduler-v5', 'Sirius-native-cpu-dispatch-admission-v6'),
]:
    new_text[name] = replace((source_root / name).read_text(), a, b)
patch = ''
after = dict(before)
modified = {}
for name, new in new_text.items():
    old = (source_root / name).read_text()
    patch += ''.join(difflib.unified_diff(old.splitlines(True), new.splitlines(True), n=0,
                                         fromfile='a/' + name, tofile='b/' + name))
    after[name] = pin(new.encode())
    modified[name] = {'basis': before[name], 'candidate': after[name]}
(provider / 'provider.patch').write_text(patch)
contract = {'scope': 'Complete direct-dispatch pure admission before LLVM. All381D functions/threeobjects unchanged. Admitted jobs use original workgroups and unchanged try_packet under existing worker FP guard; no JIT generated. All legitimate refusals use original eligible intergroup key/JIT before queue/output. Pending variant validity survives admitted launches. Coarse-ineligible/indirect ordering, original inputs/arithmetic/counts/query/pool/wait contracts retained. Unexpected worker false after immutable full admission is an internal invariant failure, explicitly aborting before NULL variant dereference; not a supported numerical/admission refusal.',
    'basis_source_identity': original['full_candidate_source_identity'],
    'modified_files': modified, 'full_candidate_source_identity': canonical(after),
    'native_objects': objects, 'patch': pin(patch.encode()),
    'source_inputs': compiled['source_inputs'],
    'admission_compile_owner': pin((work / 'build/compile/owner.json').read_bytes()),
    'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp',
    'limits': 'Source/build preparation only. Does not establish C3 timeout cause or predict acceptance. Provider/C8/soleABBA/all8/frame/science/release remain unproved.'}
(provider / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
print(json.dumps({'patch_written': True, 'shared_source_modified': False,
                  'changed_provider_files': len(modified), 'original_native_objects_reused': 3,
                  'new_pure_admission_objects': 1}))
