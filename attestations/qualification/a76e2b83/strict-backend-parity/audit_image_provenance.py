"""Verify the original three #100 reconstruction bundles from immutable Git blobs."""

import hashlib
import json
from pathlib import Path
import subprocess


REPOSITORY = 'sunyaisblank/Sirius'
BUNDLE_ROOT = 'attestations/forensic-images/b739f92'
PAIRS = (
    ('kerr-png', 'forensic-kerr-calibration', 'png'),
    ('kerr-exr', 'forensic-kerr-exr', 'exr'),
    ('alcubierre-png', 'forensic-alcubierre-completion', 'png'),
)


def read_blob(revision, path):
    return subprocess.check_output(['git', 'cat-file', 'blob', f'{revision}:{path}'])


def hash_blob(revision, path):
    sha = hashlib.sha256()
    size = 0
    with subprocess.Popen(['git', 'cat-file', 'blob', f'{revision}:{path}'],
                          stdout=subprocess.PIPE) as git:
        for chunk in iter(lambda: git.stdout.read(1 << 16), b''):
            size += len(chunk)
            sha.update(chunk)
        if git.wait() != 0:
            raise RuntimeError(f'Could not read {revision}:{path}')
    return {'path': path, 'bytes': size, 'sha256': sha.hexdigest()}


def verify_pair(name, comparison_folder, extension, remote_tags):
    tag = f'refs/tags/evidence/issue-100-{name}-reconstruction-pair'
    revision = subprocess.check_output(['git', 'rev-parse', tag], text=True).strip()
    if remote_tags.get(tag) != revision:
        raise RuntimeError(f'Remote preservation differs for {tag}')
    bundle = f'{BUNDLE_ROOT}/{name}-reconstruction-pair'
    manifest = json.loads(read_blob(revision, f'{bundle}/manifest.json'))
    report_path = f'{bundle}/files/out/{comparison_folder}/pair-result.json'
    report_bytes = read_blob(revision, report_path)
    report = json.loads(report_bytes)
    metrics = report['comparison']
    if report['evidence_kind'] != 'NEW_RECONSTRUCTION' or report['historical_tape_identity_claimed']:
        raise RuntimeError(f'Wrong evidence classification for {name}')
    if not (metrics['equal'] and metrics['encoded_files_equal'] and
            metrics['mismatched_pixels'] == metrics['mismatched_samples'] == 0 and
            metrics['compared_pixels'] == metrics['width'] * metrics['height'] and
            metrics['compared_samples'] == metrics['compared_pixels'] * metrics['channels']):
        raise RuntimeError(f'Invalid comparison for {name}')
    payloads = {}
    for producer in ('legacy', 'ported'):
        relative = f'renders/forensic-reconstruction-b739f92/{producer}/{name}.{extension}'
        archived_path = f'files/{relative}'
        blob = hash_blob(revision, f'{bundle}/{archived_path}')
        original, = [item for item in report['inputs'] if item['path'].endswith('/' + relative)]
        archived = manifest['files'][archived_path]
        if any(record[key] != blob[key] for record in (original, archived) for key in ('bytes', 'sha256')):
            raise RuntimeError(f'Recorded image differs for {name} {producer}')
        payloads[producer] = blob
    if payloads['legacy']['sha256'] != payloads['ported']['sha256']:
        raise RuntimeError(f'Encoded images differ for {name}')
    return {
        'tag': tag, 'revision': revision, 'source_export_revision': report['source_export_revision'],
        'legacy_runtime_origin': report['legacy_runtime_origin'], 'effective_spp': report['effective_spp'],
        'comparison': metrics, 'comparison_report': {
            'path': report_path, 'bytes': len(report_bytes),
            'sha256': hashlib.sha256(report_bytes).hexdigest(),
        },
        'images': payloads,
    }


if __name__ == '__main__':
    refs = subprocess.check_output(
        ['git', 'ls-remote', '--tags', 'origin', 'refs/tags/evidence/issue-100-*'], text=True)
    tags = {tag: sha for sha, tag in (line.split() for line in refs.splitlines())}
    images = {name: verify_pair(name, folder, extension, tags) for name, folder, extension in PAIRS}
    out = {
        'issue': f'https://github.com/{REPOSITORY}/issues/100',
        'classification': 'NEW_RECONSTRUCTION',
        'original_tape_identity_established': False,
        'verified_pair_payloads': images,
    }
    path = Path('out/strict-backend-science/image-recovery-audit.json')
    path.write_text(json.dumps(out, indent=2) + '\n')
    for name, pair in images.items():
        print(name, pair['comparison']['compared_pixels'], pair['images']['legacy']['sha256'])
