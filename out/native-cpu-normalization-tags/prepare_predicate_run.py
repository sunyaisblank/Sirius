"""Seal one predicate-only run after explicit admission of the existing link."""
from pathlib import Path
import hashlib,json,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;P=W/'predicate';J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
assert not subprocess.check_output(['git','status','--porcelain']).strip();actions=J(P/'actions.json');assert 'run' not in actions and not (P/'run').exists()
b=J(P/'link-result.json');assert b['passed'] and b['new_compilations']==0 and b['new_links']==1 and b['original_link_owner_passed'] is False
assert b['original_link_owner']==pin(P/'link/owner.json') and b['verification_reader']==pin(W/'admit_predicate_link.py')
inputs=dict(J(P/'prepared.json')['inputs']);inputs.update(b['runtime_libraries']);inputs[b['executable']['path']]={k:b['executable'][k] for k in ('bytes','sha256')}
for p in (P/'prepared.json',P/'link-result.json',P/'link/owner.json',P/'link/stderr.log',W/'admit_predicate_link.py',Path(__file__).resolve()):inputs[str(p.relative_to(R))]=pin(p)
for n,v in inputs.items():assert pin(R/n)==v,n
a={'revision':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'kind':'diagnostic','expected_exit':0,'software_provider':False,'limits':{'seconds':300,'rss_kib':4194304,'reserve_kib':2097152},'frozen_inputs':sorted(inputs),'argv':['/usr/bin/python3','-B',str(W/'predicate_task.py'),'run'],'scope':'One predicate-only workload on the existing verified ELF; original failed postlink owner remains failed. No compiler/linker or arithmetic workload repeats. Compare exact shared ISPC tag predicate to independent std::nextafter reference across49tuples/broadcast/mixed/masks/columns, preserving inactive tags and input words. No speed/science/release claim.'}
actions['run']=a;(P/'actions.json').write_text(json.dumps(actions,indent=2)+'\n');print(json.dumps({'prepared':'predicate-only-run','recompiled':False,'relinked':False}))
