"""Seal one original-recipe raw3 compile after the direct-bank ABI comparison."""
from pathlib import Path
import hashlib,json,os,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;D=R/'bin/linux-gcc/native-cpu-normalization-tags/objects'
J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
assert not os.environ.get('ISPC_ARGS') and not subprocess.check_output(['git','status','--porcelain']).strip()
assert not (W/'prepared-ispc-dense-inputs.json').exists()
assert D.is_dir() and {p.name for p in D.iterdir()}=={'rows_adapter.o'}
cpp=J(W/'cpp-build-result.json');assert cpp['passed'] and cpp['new_relocatable_links']==1
for phase in ('build','run'):
 owner=J(W/'cpp'/phase/'owner.json');assert owner['passed'] and owner['status']=='terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
 assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
probe=J(W/'cpp/run/stdout.log');assert probe['passed'] and probe['calls']==370440 and probe['compared_words']==44452800 and probe['tag_words']==8890560
assert J(W/'cpp-code-inspection.json')['passed'] and cpp['adapter_input_code_association']['passed']
subprocess.run(['/usr/bin/python3','-B',str(W/'derive_ispc.py')],check=True)
derived=J(W/'ispc-derivation.json');assert derived['instructions']==9380 and all(v['whole_source_inverse_restores_AI'] for v in derived['stages'].values())
old_path=R/'out/native-cpu-direct-compound/ispc-dense-compile-result.json';old=J(old_path)
assert old['all_objects_compiled'] and old['all_prepared_source_bytes_unchanged'] and old['new_compilations']==3
flags=['-O2','--target=avx2-i32x8','--cpu=core-avx2','--arch=x86-64','--addressing=32','--PIC','--opt=disable-fma','--math-lib=system']
inputs={**derived['source_inputs'],**J(W/'prepared-cpp.json')['inputs']};objects=[]
for row in old['objects']:
 slug=row['stage'];entry=derived['stages'][slug]
 assert row['exit']==0 and row['argv_executed'] and row['argv'][2:-3]==flags and row['argv'][-3:-1]==['--emit-obj','-o']
 assert row['argv'][1]==str(R/entry['basis_path'])
 args=list(row['argv']);args[1]=str(R/entry['derived_path']);args[-1]=str(D/(slug+'.o'))
 inverse=list(args);inverse[1],inverse[-1]=row['argv'][1],row['argv'][-1];assert inverse==row['argv']
 inputs[entry['derived_path']]=entry['derived']
 objects.append({'stage':slug,'output':str((D/(slug+'.o')).relative_to(R)),'argv':args,'basis_AG_argv':row['argv'],'argv_inverse_exact':True})
predicate_argv=list(objects[0]['argv']);predicate_argv[1]=str(W/'predicate.ispc');predicate_argv[-1]=str(D/'predicate.o')
compiler=Path(objects[0]['argv'][0]);assert pin(compiler)['sha256']=='1813a6b83442dfcc23a1a82a641f22452fedbff7b429877e2f31957957e0168e'
inputs[str(compiler.relative_to(R))]=pin(compiler)
resolution=subprocess.check_output(['/usr/bin/ldd','-r',str(compiler)],text=True,stderr=subprocess.STDOUT);assert 'undefined symbol' not in resolution and 'not found' not in resolution
for name in re.findall(r'(/[^\s()]+)',resolution):
 p=Path(name)
 if p.is_file():p=p.resolve(strict=True);inputs[str(p)]=pin(p)
for p in (old_path,W/'prepared-cpp.json',W/'cpp-build-result.json',W/'cpp-code-inspection.json',W/'cpp/run/owner.json',W/'cpp/run/stdout.log',*[W/n for n in ('derive_ispc.py','ispc-derivation.json','prepare_ispc_compile.py','build_ispc.py','inspect_objects.py','inspect_dense.py','compile/run_owned.py','compile/owned_guard.py','normalization_tag.isph','predicate.ispc')],R/'out/native-cpu-simd-rows/inspect_objects.py',Path('/usr/bin/nm'),Path('/usr/bin/objdump')):
 inputs[str(p.relative_to(R)) if p.is_relative_to(R) else str(p)]=pin(p)
for row in [cpp['merged_adapter'],*cpp['adapter_input_objects']]:
 p=R/row['path'];assert pin(p)=={k:row[k] for k in ('bytes','sha256')};inputs[row['path']]=pin(p)
for n,v in inputs.items():assert pin(R/n)==v,n
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
scope='Compile AJ raw3 and one predicate-only fixture once with original raw compiler/flags, after370440call/44452800payload-word/8890560tag-word scalar comparison. Shared pure-integer predicate exactly matches the unchanged AH domain. Both generic and uniform partial writes invalidate metadata, complete writes/direct CPP stores refresh it; full context initialization includes it. Original9380tuples and all arithmetic, public layouts, guards, precision, FP/ISA, requests and budgets remain. Reuse current already joined CPP/original I adapter. No new Sirius host/shader or C++ compilations/links in this action. Actual predicate/code/full38/C8/control5/frame still required; no speed/science/release claim.'
p={'revision':revision,'inputs':inputs,'objects':objects,'fixed_flags':flags,'compiler_runtime_resolution':resolution,'scope':scope,'merged_adapter':cpp['merged_adapter'],'adapter_input_code_association':cpp['adapter_input_code_association'],'adapter_input_objects':cpp['adapter_input_objects'],'predicate_argv':predicate_argv,'predicate_output':str((D/'predicate.o').relative_to(R))}
(W/'prepared-ispc-dense-inputs.json').write_text(json.dumps(p,indent=2)+'\n')
a={'revision':revision,'kind':'build_only','expected_exit':0,'software_provider':False,'limits':{'seconds':300,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/'build_ispc.py')],'frozen_inputs':sorted(set(inputs)|{str((W/'prepared-ispc-dense-inputs.json').relative_to(R))}),'scope':scope}
(W/'compile/actions.json').write_text(json.dumps({'ispc-dense':a},indent=2)+'\n');print(json.dumps({'new_ISPC_compilations_planned':4,'new_CPP_compilations_planned':0,'executed':False}))
