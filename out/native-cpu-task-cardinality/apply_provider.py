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
assert len(contract['modified_files']) == 4 and len(contract['native_objects']) == 8
assert set(contract['modified_files']) == {
    'src/gallium/drivers/llvmpipe/lp_cs_tpool.c', 'src/gallium/drivers/llvmpipe/lp_state_cs.c',
    'src/gallium/frontends/lavapipe/lvp_device.c', 'src/gallium/drivers/llvmpipe/lp_screen.c'}
subprocess.run(['/usr/bin/patch', '--dry-run', '--batch', '--fuzz=0', '-p1', '-i', str(P / 'provider.patch')], cwd=S, check=True)
subprocess.run(['/usr/bin/patch', '--batch', '--fuzz=0', '-p1', '-i', str(P / 'provider.patch')], cwd=S, check=True)
assert identity() == contract['full_candidate_source_identity']
review = {'passed': True, 'source_contract': pin((P / 'source-contract.json').read_bytes()),
    'patch': contract['patch'], 'full_source_identity_verified_before_and_after': True,
    'scope': 'Coordinator exact source review/application only. All8 X native object bytes and cs_exec_fn FP/worker body unchanged; one admitted task uses that same worker directly, and notifications follow positive task cardinality. Original thread/packet/query/count/field/numerical/fallback contracts unchanged. Not independent or machine-checked scheduling proof; requires original gates.'}
(P / 'review-admission.json').write_text(json.dumps(review, indent=2) + '\n')
print(json.dumps({'applied_files': 4, 'full_source_identity_verified': True, 'runtime_qualified': False}))
