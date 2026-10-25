"""Build/run one extended reference fixture with the exact already-built helper."""
from pathlib import Path
import hashlib,json,re,shlex,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/native-cpu-normalized-arithmetic/domain'
doc=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
mode,=sys.argv[1:];assert mode in ('build','run')
assert not subprocess.check_output(['git','status','--porcelain']).strip()
old=doc(W/'prepared-cpp.json');built=doc(W/'cpp-build-result.json');owner=doc(W/'cpp/build/owner.json')
assert built['passed'] and owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
helper=next(x['object'] for x in built['compilations'] if x['name']=='retained8')
inputs=dict(old['inputs']);inputs[helper['path']]=pin(R/helper['path'])
for p in [W/'prepared-cpp.json',W/'cpp-build-result.json',W/'cpp/build/owner.json',W/'cpp-code-inspection.json',W/'domain_driver.cpp',W/'prepare_domain.py',W/'build_domain.py',W/'domain/run_owned.py',W/'domain/owned_guard.py']:
 inputs[str(p.relative_to(R))]=pin(p)
if mode=='build':
 assert not B.exists() and not (W/'prepared-domain.json').exists();B.mkdir(parents=True)
 entry=next(x for x in old['compilations'] if x['name']=='driver');args=list(entry['argv']);args[args.index('-c')+1]=str(W/'domain_driver.cpp');args[args.index('-o')+1]=str(B/'driver.o')
 flags=args[:args.index('-c')];deps=subprocess.check_output([*flags,'-M','-MT','domain-fixture',str(W/'domain_driver.cpp')],text=True)
 for name in shlex.split(deps.replace('\\\n',' ').split(':',1)[1]):
  p=Path(name).resolve();inputs[str(p)]=pin(p)
 prepared={'revision':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'inputs':inputs,'compile':args,'link':[args[0],str(B/'driver.o'),str(R/helper['path']),'-o',str(B/'abi-check')],'helper':helper,'executable':str((B/'abi-check').relative_to(R))}
 (W/'prepared-domain.json').write_text(json.dumps(prepared,indent=2)+'\n');argv=['/usr/bin/python3','-B',str(W/'build_domain.py')]
else:
 prepared=doc(W/'prepared-domain.json');r=doc(W/'domain-build-result.json');o=doc(W/'domain/build/owner.json');assert r['passed'] and o['passed'] and o['source_unchanged'] and o['whole_inputs_unchanged']
 inputs.update(prepared['inputs']);inputs.update(r['runtime_libraries']);inputs[prepared['executable']]=r['executable_identity']
 for p in [W/'domain-build-result.json',W/'domain/build/owner.json']:inputs[str(p.relative_to(R))]=pin(p)
 argv=[str(R/prepared['executable'])]
for n,v in inputs.items():assert pin(R/n)==v,n
inputs[str((W/'prepared-domain.json').relative_to(R))]=pin(W/'prepared-domain.json')
actions=doc(W/'domain/actions.json') if (W/'domain/actions.json').exists() else {};assert mode not in actions
actions[mode]={'revision':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'kind':'build_only' if mode=='build' else 'diagnostic','expected_exit':0,'limits':{'seconds':300,'rss_kib':4194304,'reserve_kib':2097152},'argv':argv,'frozen_inputs':sorted(inputs),'software_provider':False,'scope':'One extended original scalar-reference fixture. Original mixed49cases plus broadcast and mixed normalized packets explicitly exercise full8-lane fast domains and existing fallback, all4operations/6masks/3destinations/inputaliases. Same compiled helper bytes, zero arithmetic/provider/ISPC recompilation. Finite whole-word agreement is not universal proof, performance or science/release qualification.'}
(W/'domain/actions.json').write_text(json.dumps(actions,indent=2)+'\n');print(json.dumps({'prepared':mode,'executed':False}))
