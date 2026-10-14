from pathlib import Path
import hashlib,json,os,shlex,subprocess,sys,time
root=Path.cwd();work=Path(__file__).resolve().parent;folder=Path(sys.argv[1]).resolve(strict=True);assert folder==work/'candidate-translation'
stock=root/'bin/linux-gcc/macos-moltenvk-unforced/stock';build=stock.parent/'candidate';assert not build.exists();build.mkdir()
accepted=json.loads((work/'observations/stock-byte-acceptance.json').read_text());assert accepted['passed'] and accepted['stock']['sha256']=='b2d7daec9928356de4ca24a724d135adc59685a330b043e4fa55f2fe2f9b5550'
owner=json.loads((work/'stock-translation/owner.json').read_text());assert owner['passed'] and owner['whole_inputs_unchanged'] and owner['source_unchanged'] and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
freeze=json.loads((work/'input-freeze.json').read_text())['files'];original={k.removeprefix('sources/cross-stock/'):v for k,v in freeze.items() if k.startswith('sources/cross-stock/')};candidate={k.removeprefix('sources/cross-patched/'):v for k,v in freeze.items() if k.startswith('sources/cross-patched/')};assert original.keys()==candidate.keys() and [k for k in original if original[k]!=candidate[k]]==['spirv_msl.cpp']
commands=json.loads((stock/'compile_commands.json').read_text());entry=next(x for x in commands if x['file'].endswith('/spirv_msl.cpp'));argv=shlex.split(entry['command']);source=work/'sources/cross-patched'
assert argv[0]=='/usr/bin/g++-14' and argv[-2]=='-c' and argv[-1]==str(work/'sources/cross-stock/spirv_msl.cpp')
argv=[x.replace(str(work/'sources/cross-stock'),str(source)) for x in argv];index=argv.index('-o');argv[index+1]=str(build/'spirv_msl.cpp.o')
libs=[stock/'libspirv-cross-glsl.a',stock/'libspirv-cross-core.a']
def identity(p):
 b=p.read_bytes();return {'path':str(p),'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
reused=[identity(p) for p in libs]
plan=[argv,['/usr/bin/ar','qc',str(build/'libspirv-cross-msl.a'),str(build/'spirv_msl.cpp.o')],['/usr/bin/ranlib',str(build/'libspirv-cross-msl.a')],['/usr/bin/g++-14','-std=c++17','-O2','-I'+str(source),str(work/'translate_initialize.cpp'),str(build/'libspirv-cross-msl.a'),*[str(p) for p in libs],'-o',str(build/'translate_initialize')],[str(build/'translate_initialize'),str(work/'captured/initialize.spv'),str(work/'observations/candidate.metal')]]
env=dict(os.environ,LC_ALL='C',GIT_CEILING_DIRECTORIES=str(root/'out'),GIT_CONFIG_NOSYSTEM='1',GIT_TERMINAL_PROMPT='0',PATH='/usr/bin:/bin')
for key in ('CC','CXX','CFLAGS','CXXFLAGS','CPPFLAGS','LDFLAGS','CPATH','CPLUS_INCLUDE_PATH','LIBRARY_PATH','LD_LIBRARY_PATH','LD_PRELOAD','LD_AUDIT','GIT_DIR','GIT_WORK_TREE','GIT_INDEX_FILE'):assert not env.get(key),key
rows=[]
for n,args in enumerate(plan):
 start=time.monotonic();r=subprocess.run(args,cwd=build,env=env,capture_output=True);(folder/f'{n:02d}.stdout').write_bytes(r.stdout);(folder/f'{n:02d}.stderr').write_bytes(r.stderr);rows.append({'argv':args,'cwd':str(build),'exit':r.returncode,'wall_seconds':time.monotonic()-start});(folder/'commands.json').write_text(json.dumps(rows,indent=2)+'\n');assert r.returncode==0,rows[-1]
assert reused==[identity(p) for p in libs]
receipt={'reused_unchanged_stock_libraries':reused,'recompiled_translation_unit':'spirv_msl.cpp','compile_flags_source':entry,'candidate_msl':identity(work/'observations/candidate.metal'),'byte_acceptance_pending':True,'provider_or_dispatch_executed':False,'qualification_claimed':False};(folder/'lineage.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
