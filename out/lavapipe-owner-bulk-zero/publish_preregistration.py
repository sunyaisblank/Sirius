"""Register the original finite gates after actual code admission, before C8."""
from pathlib import Path
import importlib.util,json,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('actual_provider_binding',W/'provider_binding.py')
binding=importlib.util.module_from_spec(spec);spec.loader.exec_module(binding);binding.admit()
assert not subprocess.check_output(['git','status','--porcelain']).strip()
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
p= json.loads((W/'providers.json').read_text())['providers']['candidate']
body=f'''Distinct bulk-zero trial preregistration: source {revision} on development/lavapipe-owner-bulk-zero-2026-10-23.

The preceding CodeGenNone6/22 and mixed C8timeout remain rejected. This trial retains original Default optimized machine generation and original IR passes. Only the exact Endpoint/DopriPhase intergroup-owner scratch-array zero initialization changes representation: original separate vector stores become a complete nonvolatile bulk clear at the same initialization point and lifetime. Every byte/lane, including padding/inactive owners, remains initialized; target ABI size/alignment, shape/hash/dispatch/IEEE32 admission, all fields and algorithms remain. No skipped initialization or undefined-payload assumption. New private UUID/cache namespace isolates objects.

Actual provider {p['library_identity']['sha256']}; one incremental build21.857724s,739sealed reusable build inputs unchanged. Same f958 mixed host and three immutable SPVs, zero new Sirius/shader/observer/probe compilation. Actual E/DP4-byte guard captures and complete optimized/code review must bind the full clear, equivalent remaining dataflow, Default machine pipeline, resolved libc publication, factory features/properties and FP restoration before C8. Ordinary cooperative Transport is outside both hash gates and reuses its unchanged source/code proof; no fresh Tcapture is claimed. Guard captures exit before retained arithmetic.

Original admission is unchanged: C8/full32original readbacks first, then exactly one A1/B1/B2/A2 with3internal repeats/slot;4whole-call gains>=5%,16inclusive-stage gains>=-5%,2provisioning-inclusiveAPI gains>=0%, independent Decimal agreement. Remaining7 original controls only afterall22PASS; original cold128x128 one-SPP Kerr frame<=300s only afterall8PASS on an actually linked same-selector frameELF. No pooling/subtraction/resampling, extra warmups/retries, cache/governor/input/scene/count/threshold change or original timeout-cause inference. Stop at a failed required gate. No public adoption or scientific/Mandatory/platform/native-window/installed-release/readiness admission.
'''
target=W/'preregistration.txt';assert not target.exists();target.write_text(body)
payload=W/'preregistration-payload.json';payload.write_text(json.dumps({'body':body})+'\n')
created=json.loads(subprocess.check_output(['gh','api','repos/sunyaisblank/Sirius/issues/44/comments','--input',str(payload)],text=True));assert created['body']==body
record=json.loads(subprocess.check_output(['gh','api','repos/sunyaisblank/Sirius/issues/comments/'+str(created['id'])],text=True));assert record['body']==body
(W/'preregistration-readback.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'source':revision,'fixed_original_protocol_registered':True,'url':record['html_url']}))
