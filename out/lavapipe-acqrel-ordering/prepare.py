"""Compose the two reviewed patches with authenticated original provider source."""
from pathlib import Path
import hashlib, io, json, subprocess, tarfile

R = Path.cwd(); W = Path(__file__).resolve().parent
S = W / 'upstream/mesa-25.2.8'
assert not S.exists()
member = 'attestations/qualification/8c2be67d/lavapipe-vector-workgroup-rejected/mesa-25.2.8.tar.xz'
upstream_commit = '43389d91ca54667d12a41ee0a33d35e596fc42bd'
blob = subprocess.check_output(['git', 'show', upstream_commit + ':' + member])
assert len(blob) == 43813260
assert hashlib.sha256(blob).hexdigest() == '097842f3e49d996868b38688db87b006f7d4541e93ce86d2f341d8b3e7be7c93'
contract = json.loads((W / 'source-contract.json').read_text())
assert hashlib.sha256((W / 'ieee32-prerequisite.patch').read_bytes()).hexdigest() == contract['prerequisite_patch_sha256']
(W / 'upstream').mkdir()
with tarfile.open(fileobj=io.BytesIO(blob), mode='r:xz') as archive:
    archive.extractall(W / 'upstream', filter='data')
def identity(p):
    return {'bytes': p.stat().st_size,
            'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
before = {str(p.relative_to(S)): identity(p) for p in S.rglob('*') if p.is_file()}
for name in ['ieee32-prerequisite.patch', 'acqrel.patch']:
    subprocess.run(['patch', '--batch', '--forward', '--no-backup-if-mismatch', '-p1',
                    '-i', str(W / name)], cwd=S, check=True)
after = {str(p.relative_to(S)): identity(p) for p in S.rglob('*') if p.is_file()}
assert before.keys() == after.keys()
paths = ['src/util/u_math.h', 'src/gallium/frontends/lavapipe/lvp_device.c',
         'src/gallium/drivers/llvmpipe/lp_state_cs.c', 'src/gallium/drivers/llvmpipe/lp_rast.c',
         'src/gallium/auxiliary/draw/draw_pt.c', 'src/gallium/auxiliary/gallivm/lp_bld_nir_soa.c']
assert {p for p in before if before[p] != after[p]} == set(paths)
p = S / paths[-1]
assert p.read_bytes() == (W / 'lp_bld_nir_soa.c').read_bytes()
assert before[paths[-1]] == contract['unchanged_release_file']
assert after[paths[-1]]['sha256'] == contract['new_file_sha256']
result = {'passed': True, 'official_tar_commit': upstream_commit, 'official_tar_member': member,
    'official_tar_sha256': hashlib.sha256(blob).hexdigest(), 'original_source_files': len(before),
    'changed_files': {p: {'before': before[p], 'after': after[p]} for p in paths},
    'prerequisite': 'Exact existing private5-file IEEE32 repair, no other changes.',
    'new_change': 'Only exact NIR_MEMORY_ACQ_REL maps to system-scope LLVM acquire-release under LP_VULKAN_IEEE32_ENABLED.',
    'Sirius_sources_shaders_tables_counts_unchanged': True,
    'build_or_execution_or_performance_claimed': False}
(W / 'source-admission.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'source_admitted': True, 'original_files': len(before), 'changed_files': len(paths)}))
