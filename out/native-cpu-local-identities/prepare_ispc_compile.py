"""Seal three raw stage compiles after the actual shared identity regression."""
from pathlib import Path
import hashlib,json,os,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;A=R/'out/native-cpu-value-fact-stores';D=R/'bin/linux-gcc/native-cpu-local-identities/objects';J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())
assert not os.environ.get('ISPC_ARGS') and not subprocess.check_output(['git','status','--porcelain']).strip()
assert not (W/'prepared-ispc-dense-inputs.json').exists() and not D.exists()
reg=J(W/'fixture/result.json');assert reg['passed'] and reg['calls']==555660 and reg['changed_operator_calls']==333396
assert J(W/'fixture/build/owner.json')['passed'] and J(W/'fixture/run/owner.json')['passed']
derived=J(W/'ispc-derivation.json');assert derived['instructions']==9380
oldpath=A/'ispc-dense-compile-result.json';old=J(oldpath);cpppath=A/'cpp-build-result.json';cpp=J(cpppath)
assert old['all_objects_compiled'] and old['new_raw_compilations']==3 and cpp['passed']
inputs={**derived['source_inputs'],**J(A/'prepared-cpp.json')['inputs']};objects=[]
flags=['-O2','--target=avx2-i32x8','--cpu=core-avx2','--arch=x86-64','--addressing=32','--PIC','--opt=disable-fma','--math-lib=system']
for row in old['objects']:
 slug=row['stage'];entry=derived['stages'][slug];assert row['exit']==0 and row['argv_executed'] and row['argv'][2:-3]==flags
 assert row['argv'][1]==str(R/entry['basis_path']) and row['argv'][-3:-1]==['--emit-obj','-o']
 argv=list(row['argv']);argv[1]=str(R/entry['derived_path']);argv[-1]=str(D/(slug+'.o'))
 inputs[entry['derived_path']]=entry['derived'];objects.append(dict(stage=slug,output=str((D/(slug+'.o')).relative_to(R)),argv=argv,basis_AM_argv=row['argv']))
compiler=Path(objects[0]['argv'][0]);assert pin(compiler)['sha256']=='1813a6b83442dfcc23a1a82a641f22452fedbff7b429877e2f31957957e0168e';inputs[str(compiler.relative_to(R))]=pin(compiler)
resolution=subprocess.check_output(['ldd','-r',str(compiler)],text=True,stderr=subprocess.STDOUT);assert 'undefined symbol' not in resolution and 'not found' not in resolution
for name in re.findall(r'(/[^\s()]+)',resolution):
 p=Path(name)
 if p.is_file():p=p.resolve();inputs[str(p)]=pin(p)
for p in [oldpath,cpppath,A/'cpp-code-inspection.json',W/'fixture/result.json',W/'fixture/run/owner.json',W/'fixture/build-result.json',W/'fixture/linked-code.json',W/'ispc-derivation.json',*[W/n for n in ['derive.py','identity_store.isph','prepare_ispc_compile.py','build_ispc.py','compile/run_owned.py','compile/owned_guard.py']]]:inputs[str(p.relative_to(R))]=pin(p)
for row in [cpp['merged_adapter'],*cpp['adapter_input_objects']]:
 assert pin(R/row['path'])=={k:row[k] for k in ('bytes','sha256')};inputs[row['path']]=pin(R/row['path'])
for name,v in inputs.items():assert pin(R/name)==v,name
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
scope='Compile three stages once with original ISPC/flags after555660scalar-reference calls including333396new Add/Subtract/Multiply identity-boundary calls. Only the same AM complete-fact identity/refusal stores can finish in ISPC before original C ABI calls. All writers/arrays/9380tuples/general arithmetic/Divide/Sqrt/fields/RHS/trials/quartics/precision/FP/ISA/count/cache/governors/scheduler remain. Actual unchanged AM CPP helper/adapter is reused. No renderer performance, frame or scientific acceptance inferred; this is a component awaiting whole-path gates.'
prepared=dict(revision=revision,inputs=inputs,objects=objects,fixed_flags=flags,compiler_runtime_resolution=resolution,scope=scope,merged_adapter=cpp['merged_adapter'],adapter_input_code_association=cpp['adapter_input_code_association'],adapter_input_objects=cpp['adapter_input_objects'])
D.mkdir(parents=True);(W/'prepared-ispc-dense-inputs.json').write_text(json.dumps(prepared,indent=2)+'\n')
action=dict(revision=revision,kind='build_only',expected_exit=0,software_provider=False,limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=['/usr/bin/python3','-B',str(W/'build_ispc.py')],frozen_inputs=sorted(set(inputs)|{str((W/'prepared-ispc-dense-inputs.json').relative_to(R))}),scope=scope)
(W/'compile/actions.json').write_text(json.dumps({'ispc-dense':action},indent=2)+'\n');print(json.dumps(dict(new_raw_compilations=3,new_CPP_compilations=0,executed=False)))
