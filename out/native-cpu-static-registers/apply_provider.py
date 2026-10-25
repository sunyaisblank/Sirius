"""Apply the inspected three-file patch, verifying complete source identities."""
from pathlib import Path
import hashlib
import json
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
P = W / 'provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def identity():
    value = {str(path.relative_to(S)): pin(path.read_bytes()) for path in S.rglob('*') if path.is_file()}
    return pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())


assert not (P / 'review-admission.json').exists()
contract = json.loads((P / 'source-contract.json').read_text())
assert identity() == contract['basis_source_identity']
assert pin((P / 'provider.patch').read_bytes()) == contract['patch']
assert len(contract['modified_files']) == 3 and len(contract['native_objects']) == 8
assert set(contract['modified_files']) == {
    'src/gallium/targets/lavapipe/meson.build',
    'src/gallium/frontends/lavapipe/lvp_device.c', 'src/gallium/drivers/llvmpipe/lp_screen.c'}
subprocess.run(['/usr/bin/patch', '--dry-run', '--batch', '--fuzz=0', '-p1', '-i', str(P / 'provider.patch')], cwd=S, check=True)
subprocess.run(['/usr/bin/patch', '--batch', '--fuzz=0', '-p1', '-i', str(P / 'provider.patch')], cwd=S, check=True)
assert identity() == contract['full_candidate_source_identity']
review = {'passed': True, 'source_contract': pin((P / 'source-contract.json').read_bytes()),
    'patch': contract['patch'], 'full_source_identity_verified_before_and_after': True,
    'scope': 'Coordinator source review and exact three-file provider patch application, not independent or machine-checked correctness. New raw3 propagate exact family-guarded register counts as compile-time uniform metadata; full-source inverse restores Z and allordered9380tuples. Arithmetic bodies, private scratch result ABI, generic fallback, per-row numerical refusal, admission, FP scope/scheduling, full output/publication/rollback retained. Fresh full38/C8 and original frame required. Prior Z broader gates remain scoped predecessor evidence; no new full science/release qualification.'}
(P / 'review-admission.json').write_text(json.dumps(review, indent=2) + '\n')
print(json.dumps({'applied_files': 3, 'full_source_identity_verified': True, 'runtime_qualified': False}))
