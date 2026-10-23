"""Register one original numerical cohort after actual Transport code admission."""
from pathlib import Path
import importlib.util,json,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('actual_provider_binding',W/'provider_binding.py')
binding=importlib.util.module_from_spec(spec);spec.loader.exec_module(binding);binding.admit()
assert not subprocess.check_output(['git','status','--porcelain']).strip()
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
p=json.loads((W/'providers.json').read_text())['providers']['candidate']
body=f'''Distinct complete Transport bulk-zero preregistration: source {revision}, development/lavapipe-transport-bulk-zero-2026-10-23.

Previous mixed cooperative-T/bulk-E/DP passed22/22 and all8 but its original frame stopped300.235003s at2827/16384 without completedXML. No frame retry. This new candidate adds only the exact complete row-T hash/type to the existing complete-zero source contract: all73440bytes, allocation, original insertion point/lifetime, masks, all46inputs/all40fields/7RHS/25coeff, original Default optimizer, modes/refusals/rollback and ordinary fallback retained. Original rejected263 row-T had expanded2295vector-zero stores; this representation was never applied to it. No result is extrapolated.

Actual provider {p['library_identity']['sha256']}; one incremental build3.098441s,858sealed reused inputs unchanged. One original4B Tguard capture6.166356s observed exact73440-byte resolved-libc clear and allzero return, actual Default engine/greedy, exact complete module/features/properties/FP restore. Strict remaining-IR reader refused printer names; a new independently reviewed passive bijection then established full67254-line/62908-definition equivalence with193renamed locals and no operation/Boolean rewrite. No capture or workload rerun. E/DP implementation is unchanged and composes its authenticated prior emitted-code proof; no current-provider E/DP capture is claimed. Authentic archived7f2 all-row host/archive and all3immutable SPVs are reused; no new Sirius/shader/probe compilation. Dual-T host selection was rejected source-only because it breaks the original7kernel-load assertion and has no proved crossover; no hidden substitution or gate-tuned route.

Admission remains original: C8/all32complete readbacks first; exactly one A1/B1/B2/A2,3internal repeats/slot;4whole-call gains>=5%,16inclusive-stage gains>=-5%,2provisioning-inclusiveAPI gains>=0%; independent Decimal. Remaining7 only after22PASS; original cold128x128 one-SPP Kerr<=300s only afterall8PASS and the matching original-object frame link. Preserve7kernel loads/14buffers and all inputs/counts/assertions/cache-disable/governors/limits. No pooling/subtraction/resampling/warmup/retry, adoption or scientific/Mandatory/platform/native-window/installed-release/readiness admission. Stop at the first failed required gate.
'''
target=W/'preregistration.txt';assert not target.exists();target.write_text(body)
payload=W/'preregistration-payload.json';payload.write_text(json.dumps({'body':body})+'\n')
created=json.loads(subprocess.check_output(['gh','api','repos/sunyaisblank/Sirius/issues/44/comments','--input',str(payload)],text=True));assert created['body']==body
record=json.loads(subprocess.check_output(['gh','api','repos/sunyaisblank/Sirius/issues/comments/'+str(created['id'])],text=True));assert record['body']==body
(W/'preregistration-readback.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'source':revision,'fixed_original_protocol_registered':True,'url':record['html_url']}))
