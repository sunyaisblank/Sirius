"""Restore authenticated all-row host artifacts; reuse no old numerical results."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import zipfile

R = Path.cwd()
W = Path(__file__).resolve().parent
D = R / 'attestations/qualification/bfa095ba/native-cpu-row-ownership-finite-outcome'
commit = '9e01dfdca7798c79927827334798c72c2774eebc'


def pin(path):
    return {'bytes': path.stat().st_size,
            'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}


manifest_path = D / 'MANIFEST.json'
assert subprocess.check_output(['git', 'show', commit + ':' + str(manifest_path.relative_to(R))]) == manifest_path.read_bytes()
manifest = json.loads(manifest_path.read_text())
old = R / 'out/lavapipe-intergroup-row-owner/reused-row-products.json'
assert pin(old) == {'bytes': 3212, 'sha256': 'a333b19b302b6de86a7069f995ef4d05eeb7d00a526485d390c89cbfd0d448c1'}
products = json.loads(old.read_text())['restored_actual_products']
archives = {}
restored = []
assert not (W / 'reused-row-products.json').exists()
for name, expected in products.items():
    saved = manifest['actual_product_paths'][name]
    assert all(saved[key] == expected[key] for key in ('bytes', 'sha256', 'mode'))
    path = R / name
    if not path.exists():
        blob = manifest['content_addressed_blobs'][saved['sha256']]
        assert not blob['reused']
        archive = D / blob['archive']
        if archive not in archives:
            assert pin(archive) == manifest['archives'][archive.name]
            assert subprocess.check_output(['git', 'show', commit + ':' + str(archive.relative_to(R))]) == archive.read_bytes()
            archives[archive] = pin(archive)
        with zipfile.ZipFile(archive) as handle:
            data = handle.read(blob['member'])
        assert len(data) == saved['bytes'] and hashlib.sha256(data).hexdigest() == saved['sha256']
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as handle:
            handle.write(data)
        os.chmod(path, saved['mode'])
        restored.append(name)
    assert pin(path) == {key: saved[key] for key in ('bytes', 'sha256')}
    assert path.stat().st_mode & 0o777 == saved['mode']
matched = json.loads((R / 'out/native-cpu-row-ownership/matched-result.json').read_text())
executable = R / matched['executable']
assert pin(executable)['sha256'] == matched['executable_sha256'] == '7f2d1fa88953b954ab116233458a4cb962470f863fda9bbb3a1cad7d4a7d6a26'
binary = executable.read_bytes()
compiler = json.loads((R / 'out/native-cpu-row-ownership/compiler-result.json').read_text())
for stage, item in compiler['modules'].items():
    module = R / item['module']
    assert pin(module) == {key: item[key] for key in ('bytes', 'sha256')}
    assert binary.count(module.read_bytes()) == 1, stage
(W / 'reused-row-products.json').write_text(json.dumps({
    'preservation_commit': commit, 'prior_source': 'bfa095baf0ff9592cf648cfd49a73d5f6814c05c',
    'prior_manifest': {'path': str(manifest_path.relative_to(R)), **pin(manifest_path)},
    'restored_actual_products': products, 'newly_restored_paths': restored,
    'authenticated_archives': {str(path.relative_to(R)): value for path, value in archives.items()},
    'complete_module_spans_in_final_ELF_once': True, 'new_CPP_or_shader_compilations': 0,
    'scope': 'Actual matching all-row host bytes restored from pushed Git. All35 other routes/public source0109 remain historical authority. No archived numerical/performance/frame result transfers to this new provider.',
}, indent=2) + '\n')
print(json.dumps({'actual_existing_products': len(products), 'newly_restored': restored,
                  'new_compilations': 0}))
