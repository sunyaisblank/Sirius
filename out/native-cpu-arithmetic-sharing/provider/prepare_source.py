"""Admit only the two emitted, checked modules to the unchanged private backend."""
from pathlib import Path
import difflib
import hashlib
import json

R = Path.cwd()
W = Path(__file__).resolve().parent
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
OLD = R / 'out/lavapipe-owner-bulk-zero'

def pin(p):
    data = p.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

def canonical(value):
    data = json.dumps(value, sort_keys=True, separators=(',', ':')).encode()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

review = json.loads((OLD / 'review-admission.json').read_text())
assert review['passed']
current = {str(p.relative_to(S)): pin(p) for p in S.rglob('*') if p.is_file()}
assert len(current) == 11100
assert canonical(current) == review['full_current_source_identity']
assert json.loads((W.parent / 'build/modules/owner.json').read_text())['passed']
assert json.loads((W.parent / 'emission-result.json').read_text())['passed']
modules = json.loads((W.parent / 'compiler-result.json').read_text())['modules']
hashes = {}
for kind in ['Endpoint', 'DopriPhase']:
    row = modules[kind]
    p = R / row['module']
    assert pin(p)['sha256'] == row['sha256']
    hashes[kind] = hashlib.sha1(p.read_bytes()).hexdigest()

names = ['src/gallium/frontends/lavapipe/lvp_pipeline.c',
         'src/gallium/frontends/lavapipe/lvp_device.c',
         'src/gallium/drivers/llvmpipe/lp_screen.c']
basis = {name: (S / name).read_text() for name in names}
p = S / names[0]
text = p.read_text()
key = '   static const unsigned char owner_sha1[][SHA1_DIGEST_LENGTH] = {\n'
assert text.count(key) == 1
rows = ''.join('      {' + ','.join('0x' + digest[i:i+2]
               for i in range(0, 40, 2)) + '},\n' for digest in hashes.values())
text = text.replace(key, key + rows)
text = text.replace('These two complete-row shaders', 'These four exact complete-row shaders')
p.write_text(text)
for name, old, new in [
    (names[1], 'SIR-IWG-25.2.8-3', 'SIR-ARITH-25.2-A'),
    (names[2], 'Sirius-intergroup-owner-bulkzero-v3', 'Sirius-intergroup-cpu-arithmetic-v1')]:
    p = S / name
    text = p.read_text()
    assert text.count(old) == 1
    p.write_text(text.replace(old, new))

patch = ''.join(''.join(difflib.unified_diff(basis[name].splitlines(True),
                     (S / name).read_text().splitlines(True),
                     fromfile='a/' + name, tofile='b/' + name)) for name in names)
(W / 'provider.patch').write_text(patch)
contract = {
    'scope': 'Additional exact two native CPU arithmetic-sharing E/DP hashes only; unchanged existing IEEE32/intergroup/bulk-zero/Default machinery.',
    'basis_full_source_identity': review['full_current_source_identity'],
    'basis_review': {'path': str((OLD / 'review-admission.json').relative_to(R)),
                     **pin(OLD / 'review-admission.json')},
    'modified_files': {name: {'basis': current[name],
                       'candidate': pin(S / name)} for name in names},
    'SPIRV_SHA1': hashes,
    'program_storage_modes_lane0_guard': 'Both emitted modules passed original emitted checks; Endpoint12272B / DP12584B private storage, no shared variables/barriers/atomics, no private initializers, original LocalSize64/NoContraction/denorm/RTE/program tables. Endpoint Project body unchanged.',
    'preconditions': 'Existing exact module hash plus original compute/main/no-specialization/subgroup/one-function/no-barrier structural gates; all new modules retain the original single logical owner.',
    'unchanged': 'Original two hashes remain admitted; both Transport hashes remain excluded. No arithmetic/clear extent/lifetime/mask/ABI/LLVM or machine optimization policy changes. Native bulk-zero shapes3060/3145 vectors remain exact.',
    'cache_isolation': 'New private pipeline UUID SIR-ARITH-25.2-A and disk namespace Sirius-intergroup-cpu-arithmetic-v1.',
    'provider_patch': pin(W / 'provider.patch'),
    'admission_limit': 'Provider build plus actual4B ordinary called-code screen only, after independent source review. No numerical/performance/frame/release claim.'}
(W / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
print(json.dumps({'authenticated_basis_files': len(current), 'new_hashes': hashes,
                  'changed_files': names, 'builds_or_workloads': 0}))
