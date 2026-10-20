from pathlib import Path
import json,os,shlex,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;SOURCE=R/'out/lavapipe-vector-workgroup/upstream/mesa-25.2.8';BUILD=R/'bin/linux-gcc/lavapipe-vector-workgroup'
mode,=sys.argv[1:]
if mode=='configure':
 assert not BUILD.exists() and json.loads((R/'out/lavapipe-vector-workgroup/upstream-admission.json').read_text())['passed']
 args=[str(R/'bin/toolchains/python/bin/meson'),'setup',str(BUILD),str(SOURCE),'-Dvulkan-drivers=swrast','-Dgallium-drivers=llvmpipe','-Dplatforms=','-Dglx=disabled','-Degl=disabled','-Dgbm=disabled','-Dopengl=false','-Dgles1=disabled','-Dgles2=disabled','-Dllvm=enabled','-Dshared-llvm=enabled','-Dllvm-orcjit=false','-Dbuild-tests=false','-Dmicrosoft-clc=disabled','-Dbuildtype=release','-Db_ndebug=true','-Dprefix='+str(BUILD/'uninstalled-prefix')]
 subprocess.run(args,cwd=R,check=True)
 metadata={'argv':args,'build_environment':{k:os.environ.get(k) for k in ['CC','CXX','PATH','LLVM_CONFIG','PYTHONDONTWRITEBYTECODE']},'global_installation_performed':False,'baseline_candidate_switch':'LP_SINGLE_VECTOR_WORKGROUP default0; candidate will compile only lp_state_cs.c with value1 and use a matched source relink.'}
 (W/'configuration-command.json').write_text(json.dumps(metadata,indent=2)+'\n')
elif mode=='baseline-object':
 rows=json.loads((BUILD/'compile_commands.json').read_text());selected=[v for v in rows if Path(v['file']).name=='lp_state_cs.c'];assert len(selected)==1
 v,=selected;args=shlex.split(v['command']);output=args[args.index('-o')+1];assert output.endswith('lp_state_cs.c.o') and 'libllvmpipe.a.p/' in output
 (W/'selected-compute-compile.json').write_text(json.dumps({'entry':v,'output':output,'all_commands':len(rows)},indent=2)+'\n')
 subprocess.run(['ninja','-C',str(BUILD),'-j2',output],check=True)
elif mode=='baseline-driver':
 subprocess.run(['ninja','-C',str(BUILD),'-j2','src/gallium/targets/lavapipe/libvulkan_lvp.so'],check=True)
else:raise ValueError(mode)
