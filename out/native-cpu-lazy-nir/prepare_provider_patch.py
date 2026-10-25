"""Defer final NIR optimization of exact native shaders until ordinary fallback."""
from pathlib import Path
import difflib
import hashlib
import json

R=Path.cwd();W=Path(__file__).resolve().parent;P=W/'provider'
S=R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
U=R/'out/native-cpu-scalar-zero/provider'


def doc(path):return json.loads(path.read_text())


def pin(data):return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}


def canonical(value):return pin(json.dumps(value,sort_keys=True,separators=(',',':')).encode())


assert not (P/'provider.patch').exists()
old=doc(U/'source-contract.json')
before={str(path.relative_to(S)):pin(path.read_bytes()) for path in S.rglob('*') if path.is_file()}
assert canonical(before)==old['full_candidate_source_identity']
failed=doc(R/'out/native-cpu-scalar-zero-cohort/controls/runtime/control-3/owner.json')
assert failed['status']=='terminal' and failed['stop_reason']=='diagnostic_timeout'
capture=R/'out/native-cpu-scalar-zero-stall/capture'
assert doc(capture/'owner.json')['passed'] and doc(capture/'observation.json')['stack_snapshot']
assert not doc(capture/'observation.json')['native_calls']
stacks=(capture/'all-thread-stacks.txt').read_text()
assert all(name in stacks for name in ('nir_block_cf_tree_next','lp_build_opt_nir','llvmpipe_finalize_nir','lvp_CreateComputePipelines'))
texts={}
name='src/gallium/drivers/llvmpipe/lp_screen.c';text=(S/name).read_text()
a='''{
   lp_build_opt_nir(nir);
   return NULL;
}'''
b='''{
#if LP_VULKAN_IEEE32_ENABLED
   /* Exact SPIR-V/layout admission has already marked these buffer-only
    * native shaders. Keep the complete NIR for a possible ordinary fallback;
    * positive native dispatch never consumes its optimized instruction body.
    */
   if (nir->info.stage == MESA_SHADER_COMPUTE &&
       nir->info.cs.llvmpipe_native_cpu_stage >= 1 &&
       nir->info.cs.llvmpipe_native_cpu_stage <= 3)
      return NULL;
#endif
   lp_build_opt_nir(nir);
   return NULL;
}'''
assert text.count(a)==1;text=text.replace(a,b,1)
a='Sirius-native-cpu-scalar-zero-v19';assert text.count(a)==1
texts[name]=text.replace(a,'Sirius-native-cpu-lazy-nir-v20',1)
name='src/gallium/drivers/llvmpipe/lp_state_cs.h';text=(S/name).read_text()
a='   bool intergroup_row_owner_dispatch;\n';assert text.count(a)==1
texts[name]=text.replace(a,a+'   bool native_nir_finalized;\n',1)
name='src/gallium/drivers/llvmpipe/lp_state_cs.c';text=(S/name).read_text()
# The existing zero allocation owns this one-time state. No global options,
# shared shader-info ABI, row scheduling or native admission condition changes.
assert 'struct lp_compute_shader *shader = CALLOC_STRUCT(lp_compute_shader);' in text
a='''   struct llvmpipe_screen *screen = llvmpipe_screen(lp->pipe.screen);

   struct lp_compute_shader_variant *variant ='''
b='''   struct llvmpipe_screen *screen = llvmpipe_screen(lp->pipe.screen);

#if LP_VULKAN_IEEE32_ENABLED
   struct nir_shader *native_nir = shader->base.ir.nir;
   if (!shader->native_nir_finalized &&
       native_nir->info.stage == MESA_SHADER_COMPUTE &&
       native_nir->info.cs.llvmpipe_native_cpu_stage >= 1 &&
       native_nir->info.cs.llvmpipe_native_cpu_stage <= 3) {
      /* Reached only when an ordinary variant is actually needed. Run the
       * identical final optimization before cache serialization and before
       * generate_compute's variant-specific clone/options/SoA lowering.
       * Existing context/variant ownership serializes this mutable state.
       */
      lp_build_opt_nir(native_nir);
      shader->native_nir_finalized = true;
   }
#endif

   struct lp_compute_shader_variant *variant ='''
assert text.count(a)==1;text=text.replace(a,b,1)
assert text.index('lp_build_opt_nir(native_nir);')<text.index('lp_cs_get_ir_cache_key(variant, ir_sha1_cache_key);')<text.index('generate_compute(lp, shader, variant);')
texts[name]=text
name='src/gallium/frontends/lavapipe/lvp_device.c';text=(S/name).read_text()
a='SIR-CPU-25.2.8U';assert text.count(a)==1;texts[name]=text.replace(a,'SIR-CPU-25.2.8V',1)
patch='';after=dict(before);modified={}
for name,text in texts.items():
    original=(S/name).read_text()
    patch+=''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),n=3,
        fromfile='a/'+name,tofile='b/'+name))
    after[name]=pin(text.encode());modified[name]={'basis':before[name],'candidate':after[name]}
objects=old['native_objects'];assert len(objects)==8
for row in objects:assert pin((R/row['path']).read_bytes())=={k:row[k] for k in ('bytes','sha256')}
inputs=dict(old['source_inputs'])
for path in (U/'source-contract.json',U/'provider-build-result.json',U/'linked-native-proof.json',
             R/'out/native-cpu-scalar-zero-compatibility/result.json',capture/'owner.json',
             capture/'observation.json',capture/'all-thread-stacks.txt',W/'prepare_provider_patch.py'):
    inputs[str(path.relative_to(R))]=pin(path.read_bytes())
for name,expected in inputs.items():assert pin((R/name).read_bytes())==expected,name
contract={'basis_source_identity':old['full_candidate_source_identity'],
    'full_candidate_source_identity':canonical(after),'modified_files':modified,'patch':pin(patch.encode()),
    'native_objects':objects,'source_inputs':inputs,'new_native_compilations':0,
    'guard_shared_byte_authority':'out/native-cpu-active-prefix/native_packet.hpp',
    'scope':'V defers only final lp_build_opt_nir for the three exact-SPIR-V/layout-qualified native shaders. Complete original NIR remains; the same optimizer runs once before first ordinary variant serialization/generation. U scalar3/S raw3/pureF/Iadapter bytes, all9380instructions, Q scheduler, native admission, FP, counts, host, SPV, cache-disable policy and governors unchanged. New source/build only; original C8/ABBA22+Decimal/all8/C3/frame and full science/release remain required.',
    'source_reasoning':[
        'Marker and layout checks occur before llvmpipe_finalize_nir. Every unmarked shader retains its original eager optimizer call.',
        'CALLOC initializes the per-compute-state flag false. Only ordinary generate_variant runs the deferred optimizer; it sets true only after completion, before NIR cache serialization and variant-specific clone/options/SoA passes.',
        'The original full NIR and all native input/output admission/fallback branches remain. Intergroup eligibility still inspects that NIR; before finalization it may conservatively decline an optimization, without excluding ordinary or native supported inputs.',
        'The15s snapshot locates an early pipeline-creation delay, not the whole300s causal share. This candidate removes identified redundant native-path work; no arithmetic diagnosis, successful C3 retry or qualification is inferred.'
    ]}
P.mkdir(exist_ok=True);(P/'provider.patch').write_text(patch)
(P/'source-contract.json').write_text(json.dumps(contract,indent=2)+'\n')
assert before=={str(path.relative_to(S)):pin(path.read_bytes()) for path in S.rglob('*') if path.is_file()}
print(json.dumps({'modified_provider_files':4,'reused_native_objects':8,'new_native_compilations':0,'shared_Mesa_modified':False}))
