from pathlib import Path
import ast,hashlib,json,subprocess,copy
root=Path.cwd();work=root/'out/forensic-reconstruction-b739f92';source=work/'source'
def identity(p):
 b=p.read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def write(p,o):p.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n')
assert subprocess.check_output(['git','status','--porcelain']).strip()==b''
assert subprocess.check_output(['git','rev-parse','HEAD']).decode().strip()=='5733117ba8ba4430d723d62294f991db423fbf4c'
initial=root/'attestations/forensic-images/b739f92/reconstruction-preparation-initial'
assert (initial/'manifest.json').is_file()
base=json.loads((work/'source-manifest.json').read_text())
for name,row in base['files'].items():assert identity(source/name)=={k:row[k] for k in ('bytes','sha256')}
name='lib/glfw/CMake/glfw3Config.cmake.in';origin='d1cb3bf56d0d675e12d2b4752550b9d8b030a672';blob='4a13a88b9eb2823117789072e3237829f124d5e2'
b=subprocess.check_output(['git','show',origin+':'+name]);assert subprocess.check_output(['git','rev-parse',origin+':'+name]).decode().strip()==blob
assert len(b)==115 and hashlib.sha256(b).hexdigest()=='0d266f3f8b0ae324404dc36a8ca8af01b5bf58d15b1ce806e9e318e895e9b755'
p=source/name;assert not p.exists();p.write_bytes(b)
overlay={'evidence_kind':'NEW_RECONSTRUCTION','source_export_revision':base['revision'],'build_scaffold_additions_only':{name:{'git_commit':origin,'git_blob':blob,**identity(p),'upstream_exact':'https://raw.githubusercontent.com/glfw/glfw/7b6aead9fb88b3623e3b3725ebb42670cbe4c579/CMake/glfw3Config.cmake.in'}},'runtime_source_or_asset_changed':False,'original_configuration_failure_preserved':'evidence/issue-100-reconstruction-preparation-initial','optional_uninstall_template_not_restored':'GLFW_INSTALL is forced OFF by unchanged historical cmake/legacy_build.cmake:98'}
write(work/'scaffold-overlay.json',overlay)
f=work/'phase.py';s=f.read_text();s=s.replace("build = root / 'bin/linux-gcc/forensic-reconstruction-b739f92'", "build = root / 'bin/linux-gcc/forensic-reconstruction-b739f92'\nif phase.endswith('-scaffold'):\n    build = build / 'scaffold-retry'\noperation = phase.removesuffix('-scaffold')")
for old in ("if phase == 'configure-renderers':","elif phase == 'build-renderers':","elif phase == 'build-comparator':"):
 assert old in s;s=s.replace(old,old.replace('phase ==','operation =='))
s=s.replace("require_prior('configure-renderers')", "require_prior('configure-renderers-scaffold' if phase.endswith('-scaffold') else 'configure-renderers')")
s=s.replace("assert 'CMAKE_DISABLE_FIND_PACKAGE_Vulkan", "assert 'GLFW_INSTALL:BOOL=OFF' in text\n    assert 'CMAKE_DISABLE_FIND_PACKAGE_Vulkan",1)
ast.parse(s);f.write_text(s)
m=json.loads((work/'preparation.json').read_text());oldbuild=str(root/'bin/linux-gcc/forensic-reconstruction-b739f92/renderers');newbuild=str(root/'bin/linux-gcc/forensic-reconstruction-b739f92/scaffold-retry/renderers')
for old in ('configure-renderers','build-renderers'):
 m['preparation_commands'][old+'-scaffold']=[[arg.replace(oldbuild,newbuild) for arg in cmd] for cmd in m['preparation_commands'][old]]
for case in m['cases']:
 for role,cmd in case['future_render_argv_review_only'].items():case['future_render_argv_review_only'][role]=[arg.replace(oldbuild,newbuild) for arg in cmd]
m['source_scaffold_overlay']=str(work/'scaffold-overlay.json');m['status']='COMPARATOR_CONTROLS_PASSED_SCAFFOLD_RETRY_PREPARED';m['initial_preparation_preserved']='evidence/issue-100-reconstruction-preparation-initial';m['root_commands']=[['/usr/bin/python3.12',str(work/'run_action.py'),name] for name in ('configure-renderers-scaffold','build-renderers-scaffold')]
write(work/'preparation.json',m)
a=json.loads((work/'actions.json').read_text());oldpins=json.loads((work/'pins.json').read_text())['files']
newfiles=[str(source/name),str(work/'scaffold-overlay.json'),str(Path(__file__).resolve())]
for old in ('configure-renderers','build-renderers'):
 key=old+'-scaffold';a[key]=copy.deepcopy(a[old]);a[key]['argv'][2]=key;a[key]['scope']='NEW_RECONSTRUCTION '+key+' with exactly one explicitly restored GLFW build template; no render, original tape, scientific or release claim';a[key]['frozen_inputs']=sorted(set(a[old]['frozen_inputs']+newfiles))
write(work/'actions.json',a)
pins={name:identity(Path(name)) for name in list(oldpins)+newfiles}
seal=hashlib.sha256(json.dumps(pins,sort_keys=True,separators=(',',':')).encode()).hexdigest()
write(work/'pins.json',{'files':pins,'file_count':len(pins),'file_records_sha256':seal})
for name,row in base['files'].items():assert identity(source/name)=={k:row[k] for k in ('bytes','sha256')}
print(json.dumps({'scaffold_overlay':overlay,'pins':len(pins),'preparation':identity(work/'preparation.json'),'pins_identity':identity(work/'pins.json')}))
