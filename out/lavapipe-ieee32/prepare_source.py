"""Apply a scoped IEEE32 implementation to authenticated private Mesa source."""
from pathlib import Path
import difflib,hashlib,io,json,shutil,subprocess,tarfile
R=Path.cwd();W=Path(__file__).resolve().parent;S=W/'upstream/mesa-25.2.8';(W/'upstream').mkdir()
member='attestations/qualification/8c2be67d/lavapipe-vector-workgroup-rejected/mesa-25.2.8.tar.xz'
blob=subprocess.check_output(['git','show','43389d91:'+member]);assert hashlib.sha256(blob).hexdigest()=='097842f3e49d996868b38688db87b006f7d4541e93ce86d2f341d8b3e7be7c93'
with tarfile.open(fileobj=io.BytesIO(blob),mode='r:xz') as archive:archive.extractall(W/'upstream',filter='data')
paths=['src/util/u_math.h','src/gallium/frontends/lavapipe/lvp_device.c','src/gallium/drivers/llvmpipe/lp_state_cs.c','src/gallium/drivers/llvmpipe/lp_rast.c','src/gallium/auxiliary/draw/draw_pt.c']
original={p:(S/p).read_text() for p in paths};changed=dict(original)
p=paths[0];needle='void\nutil_fpstate_set(unsigned fpstate);';assert changed[p].count(needle)==1
changed[p]=changed[p].replace(needle,needle+'''

/* This private Vulkan-only build has no D3D/OpenGL denormal-flush contract.
 * x86-64 SSE2 is mandatory. Mask shader exceptions, select nearest-even and
 * preserve both input/output subnormals; callers restore their saved MXCSR.
 * No RTZ/flush mode is advertised, so all admitted explicit modes agree.
 */
#if defined(LP_VULKAN_IEEE32) && LP_VULKAN_IEEE32 && DETECT_ARCH_X86_64
#define LP_VULKAN_IEEE32_ENABLED 1
static inline void
util_fpstate_set_ieee32(unsigned saved_fpstate)
{
   util_fpstate_set((saved_fpstate | 0x1f80u) & ~0xe040u);
}
#else
#define LP_VULKAN_IEEE32_ENABLED 0
#endif
''')
p=paths[1];needle='#include "util/u_thread.h"';assert changed[p].count(needle)==1;changed[p]=changed[p].replace(needle,needle+'\n#include "util/u_math.h"')
needle='.shaderDenormPreserveFloat32 = false,';assert changed[p].count(needle)==1;changed[p]=changed[p].replace(needle,'.shaderDenormPreserveFloat32 = LP_VULKAN_IEEE32_ENABLED,')
needle='lvp_queue_submit(struct vk_queue *vk_queue,';assert changed[p].count(needle)==1;changed[p]=changed[p].replace(needle,'lvp_queue_submit_work(struct vk_queue *vk_queue,')
needle='static VkResult\nlvp_queue_init(';assert changed[p].count(needle)==1
changed[p]=changed[p].replace(needle,'''/* Cover direct draw/geometry execution and all error returns on the queue
 * worker. Separate compute/raster workers establish the same state below.
 */
static VkResult
lvp_queue_submit(struct vk_queue *vk_queue, struct vk_queue_submit *submit)
{
#if LP_VULKAN_IEEE32_ENABLED
   unsigned fpstate = util_fpstate_get();
   util_fpstate_set_ieee32(fpstate);
#endif
   VkResult result = lvp_queue_submit_work(vk_queue, submit);
#if LP_VULKAN_IEEE32_ENABLED
   util_fpstate_set(fpstate);
#endif
   return result;
}

'''+needle)
needle='"Mesa " PACKAGE_VERSION MESA_GIT_SHA1\n';assert changed[p].count(needle)==1
changed[p]=changed[p].replace(needle,needle+'''#if LP_VULKAN_IEEE32_ENABLED
            " [Sirius scoped IEEE32 prototype]"
#endif
''')
p=paths[2];needle='#include "util/u_memory.h"';assert changed[p].count(needle)==1;changed[p]=changed[p].replace(needle,needle+'\n#include "util/u_math.h"')
needle='cs_exec_fn(void *init_data, int iter_idx, struct lp_cs_local_mem *lmem)\n{';assert changed[p].count(needle)==1
changed[p]=changed[p].replace(needle,needle+'''
#if LP_VULKAN_IEEE32_ENABLED
   /* Applies equally to pool workers and the zero-worker direct path. */
   unsigned fpstate = util_fpstate_get();
   util_fpstate_set_ieee32(fpstate);
#endif
''')
needle='                         &thread_data);\n}';assert changed[p].count(needle)==1
changed[p]=changed[p].replace(needle,'''                         &thread_data);
#if LP_VULKAN_IEEE32_ENABLED
   util_fpstate_set(fpstate);
#endif
}''')
p=paths[3]
if '#include "util/u_math.h"' not in changed[p]:changed[p]=changed[p].replace('#include "lp_flush.h"','#include "lp_flush.h"\n#include "util/u_math.h"') if '#include "lp_flush.h"' in changed[p] else changed[p].replace('#include "lp_rast.h"','#include "lp_rast.h"\n#include "util/u_math.h"')
assert '#include "util/u_math.h"' in changed[p]
needle='util_fpstate_set_denorms_to_zero(fpstate);';assert changed[p].count(needle)==2
changed[p]=changed[p].replace(needle,'''/* Retain ordinary Gallium behavior unless this Vulkan-only build
       * implements its IEEE32 execution environment. */
#if LP_VULKAN_IEEE32_ENABLED
      util_fpstate_set_ieee32(fpstate);
#else
      util_fpstate_set_denorms_to_zero(fpstate);
#endif''')
needle='''#ifdef _WIN32
   util_semaphore_signal(&task->exited);
#endif

   return 0;''';assert changed[p].count(needle)==1
changed[p]=changed[p].replace(needle,'''#if LP_VULKAN_IEEE32_ENABLED
   util_fpstate_set(fpstate);
#endif
'''+needle)
p=paths[4];assert '#include "util/u_math.h"' in changed[p]
needle='   util_fpstate_set_denorms_to_zero(fpstate);';assert changed[p].count(needle)==1
changed[p]=changed[p].replace(needle,'''#if LP_VULKAN_IEEE32_ENABLED
   /* Vertex/geometry execution must share the Vulkan precision contract. */
   util_fpstate_set_ieee32(fpstate);
#else
   util_fpstate_set_denorms_to_zero(fpstate);
#endif''')
patch=[];pins={}
for path in paths:
    assert original[path]!=changed[path]
    (S/path).write_text(changed[path]);destination=W/'changed'/path;destination.parent.mkdir(parents=True,exist_ok=True);destination.write_text(changed[path])
    patch.extend(difflib.unified_diff(original[path].splitlines(True),changed[path].splitlines(True),fromfile='a/'+path,tofile='b/'+path))
    pins[path]={'original_sha256':hashlib.sha256(original[path].encode()).hexdigest(),'candidate_sha256':hashlib.sha256(changed[path].encode()).hexdigest()}
(W/'candidate.patch').write_text(''.join(patch))
(W/'source-admission.json').write_text(json.dumps({'passed':True,'upstream_preservation_revision':'43389d91ca54667d12a41ee0a33d35e596fc42bd','official_tar_member':member,'official_tar_sha256':hashlib.sha256(blob).hexdigest(),'changed_files':pins,'candidate_build_flag':'LP_VULKAN_IEEE32=1','feature_condition':'Vulkan-only build and x86-64 SSE2; original paths retained otherwise.','other_float_properties_unchanged':True,'integer_sqrt_divide_and_workgroup_lowering_unchanged':True,'Sirius_public_sources_and_all35_shaders_unchanged':True,'build_or_runtime_qualification':False},indent=2)+'\n')
for name in ['run_owned.py','owned_guard.py']:shutil.copyfile(R/'out/lavapipe-native-fp32'/name,W/name)
print(json.dumps({'authenticated_source':True,'changed_files':len(pins),'no_provider_built_or_installed':True}))
