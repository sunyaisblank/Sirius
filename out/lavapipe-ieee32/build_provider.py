from pathlib import Path
import json,os,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;S=W/'upstream/mesa-25.2.8';B=R/'bin/linux-gcc/lavapipe-ieee32'
assert not B.exists() and json.loads((W/'source-admission.json').read_text())['passed']
args=[str(R/'bin/toolchains/python/bin/meson'),'setup',str(B),str(S),'-Dvulkan-drivers=swrast','-Dgallium-drivers=llvmpipe','-Dplatforms=','-Dglx=disabled','-Degl=disabled','-Dgbm=disabled','-Dopengl=false','-Dgles1=disabled','-Dgles2=disabled','-Dllvm=enabled','-Dshared-llvm=enabled','-Dllvm-orcjit=false','-Dbuild-tests=false','-Dmicrosoft-clc=disabled','-Dbuildtype=release','-Db_ndebug=true','-Dc_args=-DLP_VULKAN_IEEE32=1','-Dprefix='+str(B/'uninstalled-prefix')]
(W/'configuration-command.json').write_text(json.dumps({'argv':args,'build_environment':{k:os.environ.get(k) for k in ['CC','CXX','PATH','LLVM_CONFIG','PYTHONDONTWRITEBYTECODE']},'global_installation_performed':False,'baseline_provider':'Exact already-built unmodified Mesa25.2.8/GCC14/LLVM21 control retained at43389d91. Not rebuilt or retried.'},indent=2)+'\n')
subprocess.run(args,check=True)
subprocess.run(['ninja','-C',str(B),'-j2','src/gallium/targets/lavapipe/libvulkan_lvp.so'],check=True)
