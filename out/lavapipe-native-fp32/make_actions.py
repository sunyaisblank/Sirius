from pathlib import Path
import hashlib,json,re,struct,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/lavapipe-native-fp32'
binding=json.loads((W/'oracle-binding.json').read_text());fixture=R/binding['fixture'];raw=fixture.read_bytes();assert hashlib.sha256(raw).hexdigest()==binding['fixture_sha256']
magic,version,n=struct.unpack_from('<3I',raw);assert version==1 and n==174268 and len(raw)==12+24*n
records=list(struct.iter_unpack('<6I',raw[12:]));given=(W/'input.bin').read_bytes();expected=(W/'expected.bin').read_bytes()
assert given==struct.pack('<I',n)+b''.join(struct.pack('<3I',*r[:3]) for r in records)
assert expected==b''.join(struct.pack('<3I',*r[3:]) for r in records)
inputs={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}|{str(fixture.relative_to(R)),'bin/linux-gcc/lavapipe-native-fp32/gpu_probe','/usr/bin/g++-14','/usr/bin/gdb','/usr/bin/python3','/usr/share/vulkan/icd.d/lvp_icd.json','/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so','/usr/lib/x86_64-linux-gnu/libvulkan.so.1.3.275','/etc/ld.so.cache','attestations/software-vulkan/b5d2abf/optimized-ir-stage-code-binding/optimizer_observer.so'}
inputs|=set(json.loads((W/'retained-endpoint-arguments.json').read_text())['tools'].values())
for path in [Path('/usr/bin/gdb'),Path('/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so')]:
    s=subprocess.check_output(['/usr/bin/ldd',str(path)],text=True);inputs|={str(Path(x).resolve(strict=True)) for x in re.findall(r'(/[^\s()]+)',s) if Path(x).is_file()}
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();actions={}
for name in ['build','capture']:
    f=inputs|({'bin/linux-gcc/lavapipe-native-fp32/'+x for x in ['native.spv','native.spvasm','float-layout.json']} if name=='capture' else set())
    actions[name]={'revision':revision,'kind':'diagnostic','expected_exit':0,'limits':{'seconds':180,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/('build.py' if name=='build' else 'run_capture.py'))]+([] if name=='build' else ['{folder}']),'frozen_inputs':sorted(f),'software_provider':name=='capture','scope':('Compile five native binary32 operators into the original scalar corpus, requesting only RTE32 and SignedZeroInfNanPreserve32. Exact original fixture/expected values and public-factory ELF reused; no new provider or device work.' if name=='build' else 'One passive default-native arithmetic observation. Record original corpus agreement or counterexample and actual per-worker MXCSR plus returned properties; no overrides or unsupported denormal mode. Owner PASS means diagnostic completed, not numerical/feature/integrated/benefit/render/release acceptance.')}
(W/'actions.json').write_text(json.dumps(actions,indent=2)+'\n')
