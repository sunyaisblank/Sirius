from pathlib import Path
import hashlib,json,shlex,shutil,struct,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/lavapipe-vector-validation';B.mkdir()
OLD=R/'attestations/software-vulkan/9bffa35/exact-sqrt-cursor-prototype'
raw=(OLD/'corpus-probe.slang').read_text().replace('portable_binary32-candidate.h','portable_binary32.h')
(W/'corpus64.slang').write_text(raw)
prefix=raw[:raw.index('[shader("compute")]')]
(W/'barrier16.slang').write_text(prefix+(W/'barrier16-body.slang').read_text())
shutil.copyfile(OLD/'word_io.h',W/'word_io.h')
s=(W/'word_io.h').read_text().replace('expected.size() != 3 * std::size_t(input[0])','expected.size() != 3 * std::size_t(input[0]) + 48u');(W/'word_io.h').write_text(s)
s=(OLD/'gpu_probe.cpp').read_text().replace('"--dispatch"','"--dispatch64"')
s=s.replace('if (std::string(argv[5]) != "--dispatch64") throw std::runtime_error("unknown argument");','const std::uint32_t local_size = std::string(argv[5]) == "--dispatch64" ? 64u : std::string(argv[5]) == "--dispatch16" ? 16u : 0u;\n        if (local_size == 0u) throw std::runtime_error("unknown argument");')
s=s.replace('(input[0] + 63u) / 64u','(input[0] + local_size - 1u) / local_size')
s=s.replace('for (std::size_t i = 0; i < expected.size(); ++i) actual[i] = ~expected[i];','for (std::size_t i = 0; i < expected.size(); ++i) actual[i] = i < 3u * input[0] ? ~expected[i] : expected[i];')
(W/'gpu_probe.cpp').write_text(s)
fixture=R/'bin/linux-gcc/tests/backend/portable_binary32_reference.bin';raw=fixture.read_bytes()
magic,version,count=struct.unpack_from('<3I',raw);assert (magic,version,count)==(0x32334250,1,174268) and len(raw)==12+24*count
inputs=[count];expected=[]
for row in struct.iter_unpack('<6I',raw[12:]):inputs.extend(row[:3]);expected.extend(row[3:])
expected.extend((0x937ac501 ^ (i*0x9e3779b9))&0xffffffff for i in range(48))
(W/'input.bin').write_bytes(struct.pack('<'+'I'*len(inputs),*inputs));(W/'expected.bin').write_bytes(struct.pack('<'+'I'*len(expected),*expected))
rows=json.loads((R/'bin/linux-gcc/compile_commands.json').read_text())
commands=[]
probe,=[x for x in rows if x['file'].endswith('/vulkan/vulkan_device.cpp')]
args=shlex.split(probe['command']);args[args.index('-o')+1]=str(B/'gpu_probe.o');args[args.index('-c')+1]=str(W/'gpu_probe.cpp');commands.append({'argv':args,'cwd':probe['directory']})
links=json.loads((R/'out/lavapipe-vector-inspect/build-commands.json').read_text())[1]
args=[str(B/'gpu_probe.o') if x.endswith('/probe.o') else str(B/'gpu_probe') if x.endswith('/short_input_jit_probe') else x for x in links['argv']];commands.append({'argv':args,'cwd':links['cwd']})
observer,=[x for x in rows if x['file'].endswith('/tests/backend/retained_compute_test.cpp')]
args=shlex.split(observer['command']);args[args.index('-o')+1]=str(B/'audited-observer.o');args[args.index('-c')+1]=str(R/'out/lavapipe-vector-control-source/audited-observer.cpp');commands.append({'argv':args,'cwd':observer['directory']})
linkstr,=json.loads((R/'out/lavapipe-vector-inspect-source/current-backend-link.json').read_text());args=shlex.split(linkstr);args=args[args.index('/usr/bin/g++-14'):args.index('&&',args.index('/usr/bin/g++-14'))]
old='tests/backend/CMakeFiles/sirius_backend_tests.dir/retained_compute_test.cpp.o';assert args.count(old)==1
args[args.index(old)]=str(B/'audited-observer.o');args[args.index('-o')+1]=str(B/'audited_backend_tests');commands.append({'argv':args,'cwd':str(R/'bin/linux-gcc')})
cache=(R/'bin/linux-gcc/CMakeCache.txt').read_text();tools={}
for key in ['SIRIUS_SLANGC','SIRIUS_SPIRV_OPT','SIRIUS_SPIRV_DIS','SIRIUS_SPIRV_VAL']:
 tools[key]=next(line.split('=',1)[1] for line in cache.splitlines() if line.startswith(key+':FILEPATH='))
for mode in ['corpus64','barrier16']:
 source=W/(mode+'.slang');spv=B/(mode+'.compiler.spv');opt=B/(mode+'.spv');asm=B/(mode+'.spvasm')
 for argv in [[tools['SIRIUS_SLANGC'],str(source),'-I',str(R/'src/sirius/kernels'),'-O0','-target','spirv','-profile','spirv_1_5','-entry','main','-stage','compute','-o',str(spv)],
 [tools['SIRIUS_SPIRV_OPT'],'--target-env=vulkan1.2','--ssa-rewrite','--eliminate-dead-code-aggressive','--preserve-bindings','--preserve-interface',str(spv),'-o',str(opt)],
 [tools['SIRIUS_SPIRV_VAL'],'--target-env','vulkan1.2',str(opt)],
 [tools['SIRIUS_SPIRV_DIS'],str(opt),'-o',str(asm)]]:commands.append({'argv':argv,'cwd':str(R)})
(W/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
(W/'corpus-binding.json').write_text(json.dumps({'source_fixture':str(fixture),'sha256':hashlib.sha256(raw).hexdigest(),'records':count,'expected_derived_from_original_fixture_only':True,'extra_immutable_output_tail_words':48,'barrier16_contract':'Permutation j=(5*i+7)%16; inverse i=13*(j-7)%16. All lanes rendezvous four times, even tail lanes. Atomics must join the exact active-lane bitmask and count. Original independently expected high/low/valid triples and all48 tail words checked.'},indent=2)+'\n')
for name in ['run_owned.py','owned_guard.py']:shutil.copyfile(R/'out/lavapipe-vector-producer'/name,W/name)
inputs={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}|{str(fixture.relative_to(R))}
inputs|={'out/lavapipe-vector-control-source/audited-observer.cpp'}
for command in commands:
 for arg in command['argv']:
  p=Path(arg);p=p if p.is_absolute() else Path(command['cwd'])/p
  if p.is_file():inputs.add(str(p.relative_to(R)) if p.is_relative_to(R) else str(p.resolve()))
a={'revision':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'kind':'diagnostic','expected_exit':0,'limits':{'seconds':180,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/'build.py')],'frozen_inputs':sorted(inputs),'software_provider':False,'scope':'Compile original arithmetic corpus64, focused barrier16 memory/atomic corpus and exact current original control observer with extra observations only. Reuse existing current Sirius objects; no production changes or dispatch.'}
(W/'actions.json').write_text(json.dumps({'build':a},indent=2)+'\n');print(json.dumps({'prepared':True,'commands':len(commands),'frozen_inputs':len(inputs)}))
