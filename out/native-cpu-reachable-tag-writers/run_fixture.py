from pathlib import Path
import hashlib,json,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;F=W/'fixture';J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())
a=J(F/'build-admission.json');p=J(F/'prepared.json');assert a['passed'] and not (F/'result.json').exists()
exe=R/p['executable'];assert pin(exe)=={k:a['executable'][k] for k in ('bytes','sha256')}
for name,v in p['inputs'].items():assert pin(R/name)==v,name
output=subprocess.check_output([str(exe)],text=True);r=json.loads(output)
assert r==dict(passed=True,fixed_writer_calls=190512,compared_bank_words=22861440,compared_tag_words=4572288,old_stale_tag_arithmetic_counterexamples=12,witness_bank_words=1440,case_tuples=49,stages=3,writer_forms_per_stage=4),r
r.update(executable=a['executable'],build_admission=pin(F/'build-admission.json'),linked_code=pin(F/'linked-code-admission.json'),
 scope='Actual old/fixed writer bodies: twelve stale-tag arithmetic counterexamples versus the original scalar implementation and190512 corrected writer calls with complete bank/tag/mask/index checks. No changed math compilation or renderer workload. Original build-reader failure remains; source/whole-program/provider/performance/science qualification are separate.')
(F/'result.json').write_text(json.dumps(r,indent=2)+'\n');print(output,end='')
