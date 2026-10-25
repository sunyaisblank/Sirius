"""Seal one actual-writer regression build/run with unchanged AK arithmetic."""
from pathlib import Path
import hashlib,json,re,shlex,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;F=W/'fixture';B=R/'bin/linux-gcc/native-cpu-reachable-tag-writers/fixture'
J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())
assert sys.argv[1:] in (['build'],['run'])
assert not subprocess.check_output(['git','status','--porcelain']).strip()
mode=sys.argv[1];revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
if mode=='build':
 assert not (F/'prepared.json').exists() and not B.exists()
 old=J(R/'out/native-cpu-normalization-tags/prepared-cpp.json')
 driver,=[x for x in old['compilations'] if x['name']=='driver'];flags=driver['argv'][:driver['argv'].index('-c')]
 ispc=list(J(R/'out/native-cpu-normalization-tags/prepared-ispc-dense-inputs.json')['objects'][0]['argv'])
 ispc[1]=str(W/'writers.ispc');ispc[-1]=str(B/'writers.o')
 inputs=dict(J(R/'out/native-cpu-predicated-refinement/prepared-cpp.json')['inputs'])
 dep=subprocess.check_output([*flags,'-M','-MT','writer-regression',str(W/'driver.cpp')],text=True)
 for name in shlex.split(dep.replace('\\\n',' ').split(':',1)[1]):
  p=Path(name).resolve(strict=True);inputs[str(p)]=pin(p)
 derive=J(W/'ispc-derivation.json')
 for row in derive['stages'].values():
  for key in ['basis_path','derived_path']:inputs[row[key]]=pin(R/row[key])
 math=J(R/'out/native-cpu-predicated-refinement/cpp-build-result.json')
 helper,=[x['object'] for x in math['compilations'] if x['name']=='retained8']
 assert pin(R/helper['path'])=={k:helper[k] for k in ('bytes','sha256')};inputs[helper['path']]=pin(R/helper['path'])
 for p in [*W.glob('*.py'),W/'driver.cpp',W/'cases.inc',W/'writers.ispc',W/'writer-fixture-source.json',W/'ispc-derivation.json',F/'run_owned.py',F/'owned_guard.py',R/'out/native-cpu-normalization-tags/normalization_tag.isph',R/'out/native-cpu-field-compatibility/elf_proof.py',R/'out/native-cpu-predicated-refinement/cpp-build-result.json']:
  inputs[str(p.relative_to(R))]=pin(p)
 for tool in [Path(ispc[0]),Path(flags[0])]:
  inputs[str(tool.resolve())]=pin(tool.resolve())
  resolution=subprocess.check_output(['ldd','-r',str(tool)],text=True,stderr=subprocess.STDOUT)
  assert 'undefined symbol' not in resolution and 'not found' not in resolution
  for name in re.findall(r'(/[^\s()]+)',resolution):
   p=Path(name)
   if p.is_file():p=p.resolve();inputs[str(p)]=pin(p)
 commands=[ispc,[*flags,'-c',str(W/'driver.cpp'),'-o',str(B/'driver.o')],[flags[0],str(B/'driver.o'),str(B/'writers.o'),str(R/helper['path']),'-o',str(B/'writer-check')]]
 p=dict(revision=revision,inputs=inputs,commands=commands,helper=helper,executable=str((B/'writer-check').relative_to(R)),
  scope='Compile one fixture with exact alpha-renamed24 old/fixed writer bodies and original complete context declarations, plus one new driver; reuse exact AK C++ math object. Check all three stages/four writer forms/masks/indices/full payload+tag words against independent std::nextafter predicate and stale-tag arithmetic witness against original scalar. No shader/provider/frame/science/benefit execution or arithmetic recompile.')
 (F/'prepared.json').write_text(json.dumps(p,indent=2)+'\n');inputs=set(inputs)|{str((F/'prepared.json').relative_to(R))}
 argv=['/usr/bin/python3','-B',str(W/'fixture_task.py'),'build']
else:
 p=J(F/'prepared.json');built=J(F/'build-result.json');owner=J(F/'build/owner.json')
 assert built['passed'] and owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
 assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
 assert not (F/'run').exists()
 inputs=set(p['inputs'])|set(built['runtime_libraries'])|{p['executable'],str((F/'prepared.json').relative_to(R)),str((F/'build-result.json').relative_to(R)),str((F/'build/owner.json').relative_to(R))}
 argv=['/usr/bin/python3','-B',str(W/'fixture_task.py'),'run']
for name in inputs:assert (R/name).is_file(),name
actions=J(F/'actions.json') if (F/'actions.json').exists() else {};assert mode not in actions
actions[mode]=dict(revision=revision,kind='build_only' if mode=='build' else 'diagnostic',expected_exit=0,software_provider=False,
 limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=argv,frozen_inputs=sorted(inputs),scope=p['scope'])
(F/'actions.json').write_text(json.dumps(actions,indent=2)+'\n');print(json.dumps(dict(prepared=mode,executed=False)))
