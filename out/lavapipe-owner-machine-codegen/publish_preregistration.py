"""Publish the fixed finite protocol after actual code admission, before C8."""
from pathlib import Path
import importlib.util, json, subprocess
R = Path.cwd(); W = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('actual_provider_binding', W / 'provider_binding.py')
binding = importlib.util.module_from_spec(spec); spec.loader.exec_module(binding); binding.admit()
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
provider = json.loads((W / 'providers.json').read_text())['providers']['candidate']
body = f'''Distinct scoped machine-codegen trial preregistration: source {revision} on development/lavapipe-owner-machine-codegen-2026-10-23.

The preceding intergroup18/22 rejection and mixed-shader C8 timeout remain rejected. The12s partial debugger sample did not establish the cause or phase of the original300s timeout. This trial changes only the engine-owned MCJIT machine-code optimization level, after the original optimized LLVM IR, for exact admitted Endpoint/DopriPhase intergroup-owner variants. Actual original optimized IR bytes must agree with the prior captures. Cooperative Transport and other ordinary variants retain their original machine pipeline. Exact LLVM21/MCJIT/IEEE32/source/shape gates, private UUID/cache namespace, original FP settings, all fields and consumer semantics remain. No global optimization flag, IR optnone/noinline/fast-math or feature mask change.

Actual staged provider {provider['library_identity']['sha256']}; one incremental build passed21.574s,739sealed reusable inputs unchanged. Same mixed host ELF f958a60563b64c79d3f2fb6da2815bde0effdbbe49f121034cfd212ab17b27ff and three immutable SPVs; zero new Sirius, shader, observer or probe compilation. Fresh ordinary guarded captures bind the allocator factories to each actual MCJIT emission, actual code/module bytes, factory features/properties and callback FP state. These4-byte guard captures exit before retained arithmetic and do not admit correctness or speed.

Fixed original admission: C8 first; after PASS exactly one A1/B1/B2/A2 with3internal repeats/slot,4whole-call>=5%,16inclusive-stage>=-5%,2provisioning-inclusiveAPI>=0%, independent Decimal agreement and all32full readbacks. Original remaining7 only after22PASS; original cold128x128 one-SPP Kerr frame<=300s only afterall8PASS, using an actually linked same-selector frameELF. No pooling, subtraction, resampling, cache/governor/input/scene/count/threshold change or unchanged failed retry. Stop at the first failed required gate. No compile-completion, speed, frame, public adoption, full scientific/Mandatory/platform/native-window/installed-release or readiness inference.
'''
target = W / 'preregistration.txt'; assert not target.exists(); target.write_text(body)
payload = W / 'preregistration-payload.json'; payload.write_text(json.dumps({'body': body}) + '\n')
created = subprocess.check_output(['gh', 'api', 'repos/sunyaisblank/Sirius/issues/44/comments', '--input', str(payload)], text=True)
record = json.loads(created); assert record['body'] == body
readback = json.loads(subprocess.check_output(['gh', 'api', 'repos/sunyaisblank/Sirius/issues/comments/' + str(record['id'])], text=True))
assert readback['body'] == body and readback['html_url'] == record['html_url']
(W / 'preregistration-readback.json').write_text(json.dumps(readback, indent=2) + '\n')
print(json.dumps({'source': revision, 'fixed_original_protocol_registered': True, 'url': readback['html_url']}))
