"""Apply the inspected five-file patch, verifying complete source identities."""
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
    'scope': 'Coordinator source review and exact patch application only, not independent or machine-checked correctness. Replaces only three raw SIMD objects after routing original Add/Multiply return expressions through the original private scratch writer. Exact source inverse restores X; operands captured before writes, same arithmetic bodies and Y FP/scheduling scope retained. U scalar, Q scheduling and V deferred NIR remain. Same original query64, pureF/B, D fallback, IEEE32 save/restore, zero-active no-op and JIT recovery; no host/SPV/count/precision/cache/governor changes. Requires original C8/full32 and later unchanged acceptance.'}
(P / 'review-admission.json').write_text(json.dumps(review, indent=2) + '\n')
print(json.dumps({'applied_files': 3, 'full_source_identity_verified': True, 'runtime_qualified': False}))
