from pathlib import Path
import json,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;S=R/'out/lavapipe-vector-workgroup/upstream/mesa-25.2.8'
rev=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
files=[str(p) for p in sorted(S.rglob('*')) if p.is_file()]
tools=['/usr/bin/gcc-14','/usr/bin/g++-14','/usr/bin/ninja','/usr/bin/cmake','/usr/bin/ar','/usr/bin/ld','/usr/lib/llvm-21/bin/llvm-config','/usr/lib/llvm-21/lib/libLLVM.so.21.1',str(R/'bin/toolchains/python/bin/python3'),str(R/'bin/toolchains/python/bin/meson')]
env={'CC':'/usr/bin/gcc-14','CXX':'/usr/bin/g++-14','PATH':str(R/'bin/toolchains/python/bin')+':/usr/lib/llvm-21/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin','LLVM_CONFIG':'/usr/lib/llvm-21/bin/llvm-config','PYTHONDONTWRITEBYTECODE':'1'}
actions={mode:{'revision':rev,'kind':'build','scope':'Private authenticated Mesa25.2.8 '+mode+' with existing GCC14/LLVM21 tools; matched baseline control lowering, no Sirius source/math/shader/work/caches/governor changes or provider/scientific/release qualification. No installation.','expected_exit':0,'software_provider':False,'build_environment':env,'limits':{'seconds':1800,'rss_kib':4194304,'reserve_kib':2097152},'argv':['python3','-B',str(W/'provider_step.py'),mode],'frozen_inputs':[str(W/'provider_step.py'),*tools,*files]} for mode in ['configure','baseline-object','baseline-driver']}
(W/'actions.json').write_text(json.dumps(actions,indent=2)+'\n')
print(json.dumps({'source_revision':rev,'provider_source_files':len(files),'existing_tools':len(tools),'modes':list(actions)}))
