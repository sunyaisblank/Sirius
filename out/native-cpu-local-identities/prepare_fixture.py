"""Seal one ISPC identity-boundary regression using unchanged AM math/driver code."""
from pathlib import Path
import hashlib,json,re,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;A=R/'out/native-cpu-value-fact-stores';B=R/'bin/linux-gcc/native-cpu-local-identities/fixture';F=W/'fixture';J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())
assert sys.argv[1:] in ([],['run']) and not subprocess.check_output(['git','status','--porcelain']).strip()
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();F.mkdir(exist_ok=True)
if not sys.argv[1:]:
 assert not (F/'prepared.json').exists() and not B.exists()
 cpp=J(A/'cpp-build-result.json');assert cpp['passed'] and J(A/'cpp-code-inspection.json')['passed']
 failed=R/'out/native-cpu-value-fact-stores-cohort/cold-frame';owner=J(failed/'frame/owner.json')
 assert owner['status']=='terminal' and not owner['passed'] and owner['stop_reason']=='diagnostic_timeout' and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
 manifest=R/'attestations/qualification'/owner['source']['revision'][:8]/'native-cpu-value-fact-stores-frame-timeout/MANIFEST.json'
 assert subprocess.check_output(['git','show','@{upstream}:'+str(manifest.relative_to(R))])==manifest.read_bytes()
 inputs={**J(A/'prepared-cpp.json')['inputs'],**J(W/'ispc-derivation.json')['source_inputs']}
 driver,=[x['object'] for x in cpp['compilations'] if x['name']=='driver'];helper,=[x['object'] for x in cpp['compilations'] if x['name']=='retained8']
 mapping={'sirius_native_tagged_'+n+'8':'local_tagged_'+n+'8' for n in ('add','subtract','multiply')}
 objcopy=Path('/usr/bin/objcopy').resolve();ispc=Path(J(A/'ispc-dense-compile-result.json')['objects'][0]['argv'][0]);linker=Path('/usr/bin/x86_64-linux-gnu-g++-14').resolve()
 tools={objcopy,ispc,linker,Path('/usr/bin/ld').resolve(),Path('/usr/bin/ldd').resolve()}
 for tool in list(tools):
  if tool.read_bytes()[:4]!=b'\x7fELF':continue
  output=subprocess.check_output(['ldd','-r',str(tool)],text=True,stderr=subprocess.STDOUT);assert 'not found' not in output and 'undefined symbol' not in output
  tools|={Path(n).resolve() for n in re.findall(r'(/[^\s()]+)',output) if Path(n).is_file()}
 for p in [*tools,manifest,A/'cpp-build-result.json',A/'cpp-code-inspection.json',R/driver['path'],R/helper['path'],W/'fixture.ispc',W/'identity_store.isph',W/'derive.py',W/'ispc-derivation.json',W/'prepare_fixture.py',W/'fixture_task.py',F/'run_owned.py',F/'owned_guard.py',R/'out/native-cpu-field-compatibility/elf_proof.py',R/'out/native-cpu-accumulate-compatibility/elf_proof.py',Path('/usr/include/elf.h')]:inputs[str(p.relative_to(R)) if p.is_relative_to(R) else str(p)]=pin(p)
 for name,v in inputs.items():assert pin(R/name)==v,name
 flags=J(A/'ispc-dense-compile-result.json')['objects'][0]['argv'][2:-3]
 commands=[[str(objcopy),*[f'--redefine-sym={a}={b}' for a,b in mapping.items()],str(R/driver['path']),str(B/'driver.o')],[str(ispc),str(W/'fixture.ispc'),*flags,'--emit-obj','-o',str(B/'fixture.o')],[str(linker),'-Wl,--no-relax',str(B/'driver.o'),str(B/'fixture.o'),str(R/helper['path']),'-o',str(B/'abi-check')]]
 prepared=dict(revision=revision,inputs=inputs,driver=driver,helper=helper,mapping=mapping,commands=commands,executable=str((B/'abi-check').relative_to(R)),scope='Compile one integer-only ISPC helper fixture; reuse exact AM C++ math and scalar-reference driver. Only three driver undefined names are changed, with inverse symbol/relocation equality and all allocated bytes proved. Test-only no-relax link keeps exact opcode association. Actual common include is used before original three helpers. Original555660calls/66679200bank/13335840metadata words cover333396changed Add/Subtract/Multiply calls, masks, aliases, noncanonical valid values/signed zeros/nonfinite/refusal and fresh/unknown/legacy metadata. No raw/provider/frame, new C++ compilation or performance/scientific qualification.')
 (F/'prepared.json').write_text(json.dumps(prepared,indent=2)+'\n');label='build';inputs=set(inputs)|{str((F/'prepared.json').relative_to(R))}
else:
 prepared=J(F/'prepared.json');built=J(F/'build-result.json');assert built['passed'] and J(F/'build/owner.json')['passed']
 label='run';inputs=set(prepared['inputs'])|{str((F/n).relative_to(R)) for n in ('prepared.json','build-result.json','build/owner.json','linked-code.json')}|set(built['runtime_libraries'])|{prepared['executable']}
action=dict(revision=revision,kind='build_only' if label=='build' else 'diagnostic',expected_exit=0,software_provider=False,limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=['/usr/bin/python3','-B',str(W/'fixture_task.py'),label],frozen_inputs=sorted(inputs),scope=prepared['scope'])
actions=J(F/'actions.json') if (F/'actions.json').exists() else {};assert label not in actions;actions[label]=action;(F/'actions.json').write_text(json.dumps(actions,indent=2)+'\n');print(json.dumps(dict(prepared=label,executed=False)))
