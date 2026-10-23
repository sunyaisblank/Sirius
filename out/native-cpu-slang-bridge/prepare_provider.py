"""Bind exact whole modules and typed native objects to the private provider."""
from pathlib import Path
import difflib
import hashlib
import json

R = Path.cwd()
W = Path(__file__).resolve().parent
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
P = W / 'provider'


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def canonical(value):
    data = json.dumps(value, sort_keys=True, separators=(',', ':')).encode()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


review = json.loads((R / 'out/lavapipe-owner-bulk-zero/review-admission.json').read_text())
current = {str(path.relative_to(S)): pin(path) for path in S.rglob('*') if path.is_file()}
assert len(current) == 11100 and canonical(current) == review['full_current_source_identity']
compiled = json.loads((W / 'batch-compile-result.json').read_text())
assert len(compiled['objects']) == 3 and all(item['exit'] == 0 for item in compiled['objects'])
assert json.loads((W / 'build/typed-three-stage-batch/owner.json').read_text())['passed']
for name, identity in compiled['source_inputs'].items():
    assert pin(R / name) == identity, name
for item in compiled['objects']:
    assert pin(R / item['object']['path']) == {k: item['object'][k] for k in ('bytes', 'sha256')}
batch = json.loads((W / 'batch-generation.json').read_text())
P.mkdir(exist_ok=False)
basis = {}


def replace(name, old, new):
    path = S / name
    if name not in basis:
        basis[name] = path.read_text()
    text = path.read_text()
    assert text.count(old) == 1, (name, old)
    path.write_text(text.replace(old, new))


def add(name, text):
    path = S / name
    assert not path.exists()
    basis[name] = ''
    path.write_text(text)


frontend = 'src/gallium/frontends/lavapipe/'
driver = 'src/gallium/drivers/llvmpipe/'
whitelist = '#pragma once\n/* Entire authenticated SPIR-V modules, not names or partial signatures. */\n'
for stage, slug in [('Transport', 'transport'), ('Endpoint', 'endpoint'), ('DopriPhase', 'dopri_phase')]:
    item = batch['modules'][stage]['SPIRV']
    module = R / item['module']
    assert pin(module) == {key: item[key] for key in ('bytes', 'sha256')}
    data = module.read_bytes()
    whitelist += 'static const unsigned char sirius_native_' + slug + '_spirv[] = {\n'
    whitelist += ''.join(','.join('0x%02x' % byte for byte in data[i:i+24]) + ',\n'
                         for i in range(0, len(data), 24))
    whitelist += '};\n'
add(frontend + 'sirius_native_modules.h', whitelist)
add(driver + 'sirius_native_bridge.h', (W / 'native_bridge.h').read_text())
add(driver + 'lp_native_cpu_bridge.h', (W / 'provider_packet.h').read_text())
replace('src/compiler/shader_info.h',
        '         bool llvmpipe_intergroup_row_owner;',
        '         bool llvmpipe_intergroup_row_owner;\n'
        '         /* Private exact-module native CPU bridge; zero selects ordinary JIT. */\n'
        '         uint8_t llvmpipe_native_cpu_stage;')
native_admission = '''#if LP_VULKAN_IEEE32_ENABLED
#include "sirius_native_modules.h"
#endif

/* Exact byte identity binds the three CPU objects to their complete SPIR-V
 * sources. Layout admission happens separately before descriptor lowering.
 */
static uint8_t
lvp_native_cpu_stage(const VkPipelineShaderStageCreateInfo *sinfo)
{
#if LP_VULKAN_IEEE32_ENABLED
   if (sinfo->stage != VK_SHADER_STAGE_COMPUTE_BIT || sinfo->flags ||
       sinfo->pSpecializationInfo ||
       vk_find_struct_const(sinfo->pNext, PIPELINE_SHADER_STAGE_REQUIRED_SUBGROUP_SIZE_CREATE_INFO) ||
       !sinfo->pName || strcmp(sinfo->pName, "main"))
      return 0;
   const void *code;
   size_t size;
   VK_FROM_HANDLE(vk_shader_module, module, sinfo->module);
   if (module) {
      if (module->nir)
         return 0;
      code = module->data;
      size = module->size;
   } else {
      const VkShaderModuleCreateInfo *minfo =
         vk_find_struct_const(sinfo->pNext, SHADER_MODULE_CREATE_INFO);
      if (!minfo)
         return 0;
      code = minfo->pCode;
      size = minfo->codeSize;
   }
   if (!code)
      return 0;
   if (size == sizeof(sirius_native_transport_spirv) &&
       !memcmp(code, sirius_native_transport_spirv, size))
      return 1;
   if (size == sizeof(sirius_native_endpoint_spirv) &&
       !memcmp(code, sirius_native_endpoint_spirv, size))
      return 2;
   if (size == sizeof(sirius_native_dopri_phase_spirv) &&
       !memcmp(code, sirius_native_dopri_phase_spirv, size))
      return 3;
#endif
   return 0;
}

static bool
lvp_native_cpu_layout(const struct lvp_pipeline_layout *layout)
{
   if (!layout || layout->vk.set_count != 1 || !layout->vk.set_layouts[0] ||
       layout->push_constant_size)
      return false;
   const struct lvp_descriptor_set_layout *set = get_set_layout(layout, 0);
   if (set->binding_count != 2 || set->size != 2 || set->dynamic_offset_count)
      return false;
   for (unsigned i = 0; i < 2; ++i) {
      const struct lvp_descriptor_set_binding_layout *binding = &set->binding[i];
      if (!binding->valid || binding->type != VK_DESCRIPTOR_TYPE_STORAGE_BUFFER ||
          binding->array_size != 1 || binding->stride != 1 ||
          binding->descriptor_index != i || binding->uniform_block_size)
         return false;
   }
   return true;
}

'''
replace(frontend + 'lvp_pipeline.c',
        'static VkResult\ncompile_spirv(', native_admission + 'static VkResult\ncompile_spirv(')
replace(frontend + 'lvp_pipeline.c',
        '      (*nir)->info.cs.llvmpipe_intergroup_row_owner =',
        '      (*nir)->info.cs.llvmpipe_native_cpu_stage = lvp_native_cpu_stage(sinfo);\n'
        '      (*nir)->info.cs.llvmpipe_intergroup_row_owner =')
replace(frontend + 'lvp_pipeline.c',
        '   if (nir->info.stage != MESA_SHADER_TESS_CTRL)\n',
        '   if (nir->info.stage == MESA_SHADER_COMPUTE && !lvp_native_cpu_layout(layout))\n'
        '      nir->info.cs.llvmpipe_native_cpu_stage = 0;\n'
        '   if (nir->info.stage != MESA_SHADER_TESS_CTRL)\n')
replace(frontend + 'lvp_pipeline.c',
        '      if (nir->info.stage == MESA_SHADER_COMPUTE)\n'
        '         nir->info.cs.llvmpipe_intergroup_row_owner = false;',
        '      if (nir->info.stage == MESA_SHADER_COMPUTE) {\n'
        '         nir->info.cs.llvmpipe_intergroup_row_owner = false;\n'
        '         nir->info.cs.llvmpipe_native_cpu_stage = 0;\n'
        '      }')
replace(driver + 'lp_state_cs.c', '   bool intergroup_row_owner;\n',
        '   bool intergroup_row_owner;\n   uint8_t native_cpu_stage;\n')
replace(driver + 'lp_state_cs.c', '\nenum {\n   CS_ARG_CONTEXT,',
        '\n#if LP_VULKAN_IEEE32_ENABLED\n#include "lp_native_cpu_bridge.h"\n#endif\n'
        '\nenum {\n   CS_ARG_CONTEXT,')
replace(driver + 'lp_state_cs.c',
        '   variant->jit_function(&job_info->current->jit_context,',
        '   bool native_handled = false;\n#if LP_VULKAN_IEEE32_ENABLED\n'
        '   native_handled = lp_native_cpu_packet(job_info, grid_x, grid_y, grid_z);\n'
        '#endif\n   if (!native_handled)\n'
        '      variant->jit_function(&job_info->current->jit_context,')
launch_gate = '''#if LP_VULKAN_IEEE32_ENABLED
   const nir_shader *native_nir = llvmpipe->cs->base.ir.nir;
   if (native_nir->info.stage == MESA_SHADER_COMPUTE &&
       native_nir->info.cs.llvmpipe_native_cpu_stage >= 1 &&
       native_nir->info.cs.llvmpipe_native_cpu_stage <= 3 &&
       !native_nir->info.workgroup_size_variable &&
       native_nir->info.workgroup_size[0] == 64 &&
       native_nir->info.workgroup_size[1] == 1 && native_nir->info.workgroup_size[2] == 1 &&
       !native_nir->info.shared_size && !native_nir->info.task_payload_size &&
       !info->indirect && !info->variable_shared_mem &&
       !info->last_block[0] && !info->last_block[1] && !info->last_block[2] &&
       !info->grid_base[0] && !info->grid_base[1] && !info->grid_base[2] &&
       !info->draw_count && !info->indirect_draw_count && !info->num_globals &&
       info->block[0] == 64 && info->block[1] == 1 && info->block[2] == 1 &&
       job_info.grid_size[0] >= 1 && job_info.grid_size[0] <= 65536 &&
       job_info.grid_size[1] == 1 && job_info.grid_size[2] == 1)
      job_info.native_cpu_stage = native_nir->info.cs.llvmpipe_native_cpu_stage;
#endif

'''
replace(driver + 'lp_state_cs.c',
        '   job_info.current = &llvmpipe->csctx->cs.current;\n\n   int num_tasks =',
        '   job_info.current = &llvmpipe->csctx->cs.current;\n\n' + launch_gate + '   int num_tasks =')
replace(frontend + 'lvp_device.c', 'SIR-IWG-25.2.8-3', 'SIR-CPU-25.2.8A')
replace(driver + 'lp_screen.c', 'Sirius-intergroup-owner-bulkzero-v3',
        'Sirius-native-cpu-slang-bridge-v1')
objects = ',\n    '.join("'" + str((R / item['object']['path']).resolve()) + "'"
                         for item in compiled['objects'])
replace('src/gallium/targets/lavapipe/meson.build',
        "  [ 'lavapipe_target.c' ],\n",
        "  [ 'lavapipe_target.c' ],\n  objects : files(\n    " + objects + '\n  ),\n')
patch = ''.join(''.join(difflib.unified_diff(basis[name].splitlines(True),
                      (S / name).read_text().splitlines(True),
                      fromfile=('a/' + name if basis[name] else '/dev/null'),
                      tofile='b/' + name)) for name in basis)
(P / 'provider.patch').write_text(patch)
after = {str(path.relative_to(S)): pin(path) for path in S.rglob('*') if path.is_file()}
contract = {
    'scope': 'Private native CPU bridge source preparation only. No configuration/build/runtime executed.',
    'basis_source_identity': review['full_current_source_identity'],
    'modified_files': {name: {'basis': current.get(name), 'candidate': after[name]} for name in basis},
    'full_candidate_source_identity': canonical(after),
    'source_files': len(after), 'native_objects': [item['object'] for item in compiled['objects']],
    'source_inputs': {str(path.relative_to(R)): pin(path) for path in
                      [W / 'prepare_provider.py', W / 'provider_packet.h', W / 'native_bridge.h',
                       W / 'batch-generation.json', W / 'batch-compile-result.json']},
    'module_admission': 'Complete byte comparison of all3 authenticated original native-row SPIR-Vs, main/compute/no specialization/subgroup; native marker zero on opaque NIR import. Unique new private UUID/disk namespace, serialized shader-info includes stage.',
    'descriptor_admission': 'Exactly one set, two scalar storage-buffer bindings0/1 at descriptor indices0/1, no push constants/dynamic/inline storage. constants[1] exact two-descriptor extent in four-byte units; descriptor buffer extents in bytes; reject output/table alias. Actual typed object checks complete tables/dimensions/alignment/nonalias before writes.',
    'execution': 'Inside original IEEE32 callback guard and original task lifetime/query counts. Direct base0 64x1x1/gridNx1x1 only. E/DP retain existing8group packet/tail; T ordinary1. All64 generated invocations/logicalgroup; original whole-row clearing/fields/math/refusals. Ordinary Default/greedy JIT remains available and handles every refusal before output writes.',
    'limits': 'Source/module/object association is conditional on sealed inputs. Runtime FP equivalence/libm/actual path/negative fallback/benefit/frame remain unproved. No production files or original gates altered.',
    'patch': pin(P / 'provider.patch'),
}
(P / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
print(json.dumps({'basis_files': len(current), 'candidate_files': len(after),
                  'changed_files': len(basis), 'compiler_or_workload_executed': False}))
