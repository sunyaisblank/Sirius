"""Reuse exact P products; seal one integrator compile/archive/two links."""
from pathlib import Path
import hashlib,json,re,shlex,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/native-cpu-owned-lower-increments'
doc=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert not B.exists() and not (W/'prepared-build.json').exists()
base_path=R/'out/native-cpu-upload-aware-staging/host-build-result.json';base=doc(base_path)
assert base['passed'] and base['new_cpp_compilations']==1 and base['all3modules_once_in_both_executables']
old_owner=doc(base_path.parent/'build/host-build/owner.json')
assert old_owner['passed'] and old_owner['source_unchanged'] and old_owner['whole_inputs_unchanged']
assert not old_owner['remaining_owned_processes'] and not old_owner['cleanup_errors']
proof_path=W/'host-source-proof.json';proof=doc(proof_path)
assert proof['full_inverse_restores_original_source'] and proof['all_other_source_bytes_unchanged']
assert pin(R/proof['source'])==proof['basis'] and pin(R/proof['derived_path'])==proof['derived']
assert subprocess.check_output(['git','show','0109ab46:'+proof['source']])==(R/proof['source']).read_bytes()
old_archive=R/base['backend_archive'];assert pin(old_archive)==base['backend_archive_identity']
old_member=subprocess.check_output(['/usr/bin/ar','p',str(old_archive),'retained_integrator.cpp.o'])
assert {'bytes':len(old_member),'sha256':hashlib.sha256(old_member).hexdigest()}==base['unchanged_archive_members']['retained_integrator.cpp.o']
inputs={str(p.relative_to(R)):pin(p) for p in [base_path,base_path.parent/'build/host-build/owner.json',proof_path,R/proof['source'],R/proof['derived_path'],old_archive,*W.glob('*.py')]}
for rows in base['unchanged_other_link_inputs'].values():
 for v in rows:
  p=R/v['path'];assert pin(p)=={k:v[k] for k in ('bytes','sha256')},v['path'];inputs[v['path']]=pin(p)
for v in base['restored_private_headers']:
 p=R/v['path'];assert pin(p)=={k:v[k] for k in ('bytes','sha256')};inputs[v['path']]=pin(p)
v=base['restored_observer'];assert pin(R/v['path'])=={k:v[k] for k in ('bytes','sha256')};inputs[v['path']]=pin(R/v['path'])
for v in base['module_joins'].values():
 p=R/v['module'];assert pin(p)=={k:v[k] for k in ('bytes','sha256')};inputs[v['module']]=pin(p)
obj=B/'matched/retained_integrator.cpp.o';archive=B/'matched/libsirius_backend.a';host=B/'matched/audited_backend_tests';frame=B/'frame/sirius_render_tests'
compile_args=list(base['commands'][0]['argv']);compile_args[compile_args.index('-c')+1]=str(R/proof['derived_path']);compile_args[compile_args.index('-o')+1]=str(obj)
flags=compile_args[:compile_args.index('-o')]
deps=subprocess.check_output([*flags,'-M','-MT','owned-lower',str(R/proof['derived_path'])],text=True)
for name in shlex.split(deps.replace('\\\n',' ').split(':',1)[1]):
 p=Path(name).resolve(strict=True);inputs[str(p)]=pin(p)
commands=[{'argv':compile_args,'cwd':base['commands'][0]['cwd']},{'argv':['/usr/bin/ar','rDs',str(archive),str(obj)],'cwd':str(R/'bin/linux-gcc')}]
for index,destination in [(2,host),(3,frame)]:
 args=list(base['commands'][index]['argv']);assert args.count(str(old_archive)) in (1,2)
 args=[str(archive) if x==str(old_archive) else x for x in args];args[args.index('-o')+1]=str(destination)
 commands.append({'argv':args,'cwd':base['commands'][index]['cwd']})
# Reuse tool discovery only; bind actual current bytes, including the current
# zlib already consumed by the AH helper build. Historical P receipts remain unchanged.
cpp_preparation=R/'out/native-cpu-normalized-arithmetic/prepared-cpp.json';cpp=doc(cpp_preparation)
for name,expected in cpp['inputs'].items():
 if name.startswith(('/usr/','/lib/','/opt/')):
  assert pin(Path(name))==expected,name;inputs[name]=expected
for name in ('/usr/bin/ar','/usr/bin/ranlib','/usr/bin/ld','/usr/bin/objdump','/usr/bin/nm'):
 p=Path(name).resolve();inputs[str(p)]=pin(p)
inputs[str(cpp_preparation.relative_to(R))]=pin(cpp_preparation)
for mode in ['matched','frame']:(B/mode).mkdir(parents=True,exist_ok=True)
prepared={'revision':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'commands':commands,'source_inputs':inputs,'source_proof':proof,'base_P_result':str(base_path.relative_to(R)),
 'object':str(obj.relative_to(R)),'backend_archive':str(archive.relative_to(R)),'executable':str(host.relative_to(R)),'frame_executable':str(frame.relative_to(R)),
 'scope':'One new host integrator compilation, one-member archive replacement and original audited-host/frame links. Reuse P compute/header/observer/products byte-for-byte; no native math/shader/compiler operation in this action. Only20duplicate lower position/beam increments become exact copies of the existing same-step lower PhasePacket; all other source bytes and calls remain. P spacing/upload changes remain prior unchanged evidence. New host must pass affected original DP control/C8 and original frame; no science/release claim.'}
(W/'prepared-build.json').write_text(json.dumps(prepared,indent=2)+'\n')
action={'revision':prepared['revision'],'kind':'build_only','expected_exit':0,'limits':{'seconds':300,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/'build_host.py')],'frozen_inputs':sorted(set(inputs)|{str((W/'prepared-build.json').relative_to(R))}),'software_provider':False,'scope':prepared['scope']}
(W/'build/actions.json').write_text(json.dumps({'host-build':action},indent=2)+'\n')
print(json.dumps({'new_CPP_compilations':1,'new_links':2,'frozen_inputs':len(inputs),'executed':False}))
