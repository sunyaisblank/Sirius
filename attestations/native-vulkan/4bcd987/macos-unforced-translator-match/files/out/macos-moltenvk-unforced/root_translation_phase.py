from pathlib import Path
import json,os,subprocess,sys,time
root=Path.cwd();work=Path(__file__).resolve().parent
folder=Path(sys.argv[1]).resolve(strict=True);assert folder==work/'stock-translation'
source=work/'sources/cross-stock';build=root/'bin/linux-gcc/macos-moltenvk-unforced/stock';assert not build.exists()
build.mkdir(parents=True)
env=dict(os.environ,LC_ALL='C',GIT_CEILING_DIRECTORIES=str(root/'out'),GIT_CONFIG_NOSYSTEM='1',GIT_TERMINAL_PROMPT='0',PATH='/usr/bin:/bin')
for key in ('CC','CXX','CFLAGS','CXXFLAGS','CPPFLAGS','LDFLAGS','CMAKE_TOOLCHAIN_FILE','CMAKE_PROJECT_INCLUDE','CMAKE_PROJECT_INCLUDE_BEFORE','CMAKE_PREFIX_PATH','CPATH','CPLUS_INCLUDE_PATH','LIBRARY_PATH','LD_LIBRARY_PATH','LD_PRELOAD','LD_AUDIT','GIT_DIR','GIT_WORK_TREE','GIT_INDEX_FILE'):
 assert not env.get(key),key
flags=['-DCMAKE_BUILD_TYPE=Release','-DCMAKE_EXPORT_COMPILE_COMMANDS=ON','-DCMAKE_C_COMPILER=/usr/bin/gcc-14','-DCMAKE_CXX_COMPILER=/usr/bin/g++-14','-DSPIRV_CROSS_STATIC=ON','-DSPIRV_CROSS_SHARED=OFF','-DSPIRV_CROSS_CLI=OFF','-DSPIRV_CROSS_ENABLE_TESTS=OFF','-DSPIRV_CROSS_ENABLE_GLSL=ON','-DSPIRV_CROSS_ENABLE_MSL=ON','-DSPIRV_CROSS_ENABLE_HLSL=OFF','-DSPIRV_CROSS_ENABLE_CPP=OFF','-DSPIRV_CROSS_ENABLE_REFLECT=OFF','-DSPIRV_CROSS_ENABLE_C_API=OFF','-DSPIRV_CROSS_ENABLE_UTIL=OFF']
commands=[['/usr/bin/cmake','-S',str(source),'-B',str(build),'-G','Ninja',*flags],['/usr/bin/cmake','--build',str(build),'--target','spirv-cross-msl','--parallel','1'],['/usr/bin/g++-14','-std=c++17','-O2','-I'+str(source),str(work/'translate_initialize.cpp'),str(build/'libspirv-cross-msl.a'),str(build/'libspirv-cross-glsl.a'),str(build/'libspirv-cross-core.a'),'-o',str(build/'translate_initialize')],[str(build/'translate_initialize'),str(work/'captured/initialize.spv'),str(work/'observations/stock.metal')]]
rows=[]
for index,argv in enumerate(commands):
 start=time.monotonic();r=subprocess.run(argv,cwd=source,env=env,capture_output=True)
 (folder/f'{index:02d}.stdout').write_bytes(r.stdout);(folder/f'{index:02d}.stderr').write_bytes(r.stderr)
 rows.append({'argv':argv,'cwd':str(source),'exit':r.returncode,'wall_seconds':time.monotonic()-start});(folder/'commands.json').write_text(json.dumps(rows,indent=2)+'\n')
 assert r.returncode==0,rows[-1]
print(json.dumps({'stock_translation_executed':True,'exact_byte_acceptance_pending':True,'provider_executed':False}),flush=True)
