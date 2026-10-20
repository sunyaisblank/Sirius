"""Reuse the original scalar oracle, public-factory ELF and passive observer."""
from pathlib import Path
import hashlib,io,json,subprocess,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/lavapipe-native-fp32';B.mkdir()
prefix='attestations/qualification/6b355fbf/checked-fp64-three-stage-rejected/'
def archive(name):return zipfile.ZipFile(io.BytesIO(subprocess.check_output(['git','show','a78c7e1b:'+prefix+name])))
with archive('diagnostics-and-trial.zip') as z:
    for name in ['run_owned.py','owned_guard.py','input.bin','expected.bin','gpu_probe.cpp','oracle-binding.json']:
        (W/name).write_bytes(z.read('out/portable-checked-fp64/'+name))
    for name in ['capture_jit.py','capture_features.py','read_features.py','feature-layout.json']:
        (W/name).write_bytes(z.read('out/portable-checked-fp64-jit/'+name))
    (W/'retained-endpoint-arguments.json').write_bytes(z.read('out/portable-checked-fp64/endpoint-arguments.json'))
    (W/'original-corpus.slang').write_bytes(z.read('out/portable-checked-fp64/corpus.slang'))
    (W/'original-probe-commands.json').write_bytes(z.read('out/portable-checked-fp64/commands.json'))
with archive('actual-products.zip') as z:
    raw=z.read('bin/linux-gcc/portable-checked-fp64/gpu_probe');(B/'gpu_probe').write_bytes(raw);(B/'gpu_probe').chmod(0o755)
source=(R/'src/sirius/kernels/portable_binary32.h').read_bytes();(W/'portable_binary32.h').write_bytes(source)
original=(W/'original-corpus.slang').read_text()
for name in ['Add','Subtract','Multiply','Divide']:
    original=original.replace('return PB32'+name+'(a, b);','return Native'+name+'(a, b);')
original=original.replace('return PB32Sqrt(a);','return NativeSqrt(a);')
native='''uint NativeCanonical(float value) {
    uint bits=asuint(value);
    return PB32Nan(bits) ? 0x7fc00000u : bits;
}
'''
for name,op in [('Add','OpFAdd'),('Subtract','OpFSub'),('Multiply','OpFMul'),('Divide','OpFDiv')]:
    native+='uint Native'+name+'''(uint a,uint b) {
    float x=asfloat(a), y=asfloat(b);
    float value=spirv_asm { '''+op+''' $$float %v $x $y; OpDecorate %v NoContraction; OpCopyObject $$float result %v; };
    return NativeCanonical(value);
}
'''
native+='uint NativeSqrt(uint a) { return NativeCanonical(sqrt(asfloat(a))); }\n'
(W/'corpus.slang').write_text(original.replace('#include "portable_binary32.h"','#include "portable_binary32.h"\n'+native,1))
# Record every worker used by this one dispatch. Read registers only; never
# set MXCSR, alter a property or execute a GDB inferior function call.
p=W/'capture_jit.py';s=p.read_text()
old="            if key not in self.seen:\n                self.seen.add(key)"
assert s.count(old)==1
s=s.replace(old,"            seen_key = (key, gdb.selected_thread().global_num)\n            if seen_key not in self.seen:\n                self.seen.add(seen_key)")
old="'worker_thread': gdb.selected_thread().global_num,"
assert s.count(old)==1
s=s.replace(old,old+" 'actual_mxcsr': reg('mxcsr'),")
p.write_text(s)
(W/'reuse.json').write_text(json.dumps({'source_preservation_revision':'a78c7e1b','probe_byte_exact_sha256':hashlib.sha256(raw).hexdigest(),'probe_CXX_revision':'c41bc441546421d8d34a61a7c4ee5ca5ae62ca82','underlying_public_CXX_logic':'0109ab462edade524667d20303507aa09adcbaf6','corpus_records':174268,'expected_words':522804,'original_expected_values_unchanged':True,'new_probe_CPP_compilations':0,'passive_observer_recompiled':False,'GDB_capture_delta':'Observe actual MXCSR at each worker first entry; no register/property writes.','changed_arithmetic':'Only five primitive operations use native binary32; all other original operations and the exact integer product/residual remain unchanged.','scope':'Default denormal behavior with supported RTE32 and SignedZeroInfNanPreserve32; no DenormPreserve request, feature mask, rendering or qualification.'},indent=2)+'\n')
