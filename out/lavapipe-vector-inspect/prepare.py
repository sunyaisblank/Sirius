"""Adapt the existing passive ordinary-JIT observer to the matched LLVM21 pair."""
from pathlib import Path
import hashlib,json,re,shlex,shutil,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;OLD=R/'out/lavapipe-vector-inspect-source'
B=R/'bin/linux-gcc/lavapipe-vector-inspect';B.mkdir()
for n in ['capture_jit.py','capture_features.py','read_features.py','feature-layout.json','probe-source.cpp']:
    shutil.copyfile(OLD/n,W/n)
llvm=Path('/usr/lib/llvm-21/lib/libLLVM.so.21.1').resolve(strict=True)
symbols=subprocess.check_output(['/usr/bin/nm','-D','--defined-only',str(llvm)],text=True)
def symbol(name):
    rows=[s for s in symbols.splitlines() if s.split()[-1].split('@@')[0]==name];assert len(rows)==1
    return int(rows[0].split()[0],16)
parent,passes=symbol('LLVMGetGlobalParent'),symbol('LLVMRunPasses')
disasm=subprocess.check_output(['/usr/bin/objdump','-d','--start-address='+hex(parent),'--stop-address='+hex(parent+16),str(llvm)],text=True)
assert '0x28(%rdi),%rax' in disasm
guard='f30f1efa554889e5488b47285dc3cccc'
s=(W/'capture_jit.py').read_text().replace('f30f1efa554889e5488b47285dc36690',guard)
s=s.replace('parent_api - 0x10cc070 + 0x473aad0',f'parent_api - {hex(parent)} + {hex(passes)}')
s=s.replace("'libLLVM.so.20.1'",repr(llvm.name))
(W/'capture_jit.py').write_text(s)
obs=R/'attestations/software-vulkan/b5d2abf/optimized-ir-stage-code-binding/optimizer_observer.c'
(W/'optimizer_observer.c').write_text(obs.read_text().replace('/lib/x86_64-linux-gnu/libLLVM.so.20.1',str(llvm)))
rows=json.loads((R/'bin/linux-gcc/compile_commands.json').read_text())
row,=[a for a in rows if a['file'].endswith('/vulkan/vulkan_device.cpp')]
compile=shlex.split(row['command'])
compile[compile.index('-o')+1]=str(B/'probe.o');compile[compile.index('-c')+1]=str(W/'probe-source.cpp')
linkstr,=json.loads((OLD/'current-backend-link.json').read_text())
args=shlex.split(linkstr);start=args.index('/usr/bin/g++-14');end=args.index('&&',start)
args=args[start:end];args=[x for x in args if not x.endswith('.cpp.o') and not x.endswith('/libgtest.a') and not x.endswith('/libgtest_main.a')]
args[args.index('-o')+1]=str(B/'short_input_jit_probe');args.insert(args.index('-o'),str(B/'probe.o'))
build=[{'argv':compile,'cwd':row['directory']},{'argv':args,'cwd':str(R/'bin/linux-gcc')},
 {'argv':['/usr/bin/gcc-14','-std=c11','-D_GNU_SOURCE','-Wall','-Wextra','-Werror','-fPIC','-shared','-O2',str(W/'optimizer_observer.c'),'-ldl','-o',str(B/'optimizer_observer.so')],'cwd':str(R)}]
(W/'build-commands.json').write_text(json.dumps(build,indent=2)+'\n')
source=R/'out/lavapipe-vector-workgroup/upstream/mesa-25.2.8'
for arm in ['baseline','candidate']:
    subprocess.run(['/usr/bin/python3','-B',str(source/'src/vulkan/util/vk_icd_gen.py'),'--api-version','1.4','--xml',str(source/'src/vulkan/registry/vk.xml'),'--sizeof-pointer','8','--lib-path',str(R/f'bin/linux-gcc/lavapipe-vector-staged/{arm}/libvulkan_lvp.so'),'--out',str(W/f'{arm}-icd.json')],check=True)
    assert json.loads((W/f'{arm}-icd.json').read_text())['ICD']['api_version']=='1.4.318'
(W/'adaptation.json').write_text(json.dumps({'passive_observer_reuse':str(OLD),'original_observer_C':str(obs),
 'LLVM':str(llvm),'LLVM_GetGlobalParent_offset':hex(parent),'LLVM_RunPasses_offset':hex(passes),'getter_16_byte_guard':guard,
 'getter_disassembly':disasm,'runtime_ICD_source':'Unmodified Mesa vk_icd_gen.py and official vk.xml,1.4.318; staged relink manifests1.4.305 remain unused historical metadata.',
 'no_driver_or_IR_instrumentation':True,'ordinary_pass_forwarded_unchanged':True},indent=2)+'\n')
for n in ['run_owned.py','owned_guard.py']:shutil.copyfile(R/'out/lavapipe-vector-producer'/n,W/n)
inputs={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
for cmd in build:
    for a in cmd['argv']:
        p=Path(a);p=p if p.is_absolute() else Path(cmd['cwd'])/p
        if p.is_file():inputs.add(str(p.relative_to(R)) if p.is_relative_to(R) else str(p.resolve()))
a={'revision':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'kind':'diagnostic','expected_exit':0,'limits':{'seconds':120,'rss_kib':2097152,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/'build.py')],'frozen_inputs':sorted(inputs),'software_provider':False,'scope':'Compile existing short-input ordinary-JIT observer with current strict GCC14 C++26 flags and exact existing backend link; compile passive LLVM21 pass-export observer. No retained arithmetic or timing.'}
(W/'actions.json').write_text(json.dumps({'build-observers':a},indent=2)+'\n')
print(json.dumps({'prepared':True,'build_commands':len(build),'frozen_inputs':len(inputs)}))
