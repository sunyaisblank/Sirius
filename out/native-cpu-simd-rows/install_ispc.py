"""Acquire the pinned official CPU compiler inside the existing preset cache.

Extract only compiler, license and bundled language documentation. No system
installation, tool invocation, shader generation or runtime workload occurs.
"""
from pathlib import Path
import hashlib
import json
import os
import tarfile
import time
import urllib.request

root = Path.cwd()
work = Path(__file__).resolve().parent
destination = root / 'bin/linux-gcc/native-cpu-simd-rows/compiler'
tag = 'v1.32.0'
name = 'ispc-v1.32.0-linux.tar.gz'
size = 122614155
digest = 'faf266bda2b702d68976ffd26bbd36c1ddba43424cc3c44f25b092cd1aa1f64a'
api = 'https://api.github.com/repos/ispc/ispc/releases/tags/' + tag


def request(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers={
        'User-Agent': 'Sirius-development', 'Accept': 'application/vnd.github+json'}), timeout=30)


assert not destination.exists()
assert not (work / 'compiler-acquisition.json').exists()
release = json.load(request(api))
assert release['tag_name'] == tag and not release['draft'] and not release['prerelease']
asset, = [a for a in release['assets'] if a['name'] == name]
assert asset['id'] == 637733371 and asset['size'] == size
assert asset['digest'] == 'sha256:' + digest
assert asset['browser_download_url'] == 'https://github.com/ispc/ispc/releases/download/' + tag + '/' + name
destination.mkdir(parents=True)
archive = destination / name
begin = time.monotonic()
sha = hashlib.sha256()
total = 0
with request(asset['browser_download_url']) as response, archive.open('xb') as output:
    while data := response.read(1024 * 1024):
        assert time.monotonic() - begin <= 300, 'bounded acquisition timed out'
        total += len(data)
        assert total <= size
        sha.update(data)
        output.write(data)
assert total == size and sha.hexdigest() == digest
members = {}
with tarfile.open(archive, 'r:gz') as package:
    for member in package:
        parts = Path(member.name).parts
        assert parts and not Path(member.name).is_absolute() and '..' not in parts
        assert parts[0] == 'ispc-v1.32.0-linux'
        relative = Path(*parts[1:])
        selected = relative == Path('bin/ispc') or (
            len(relative.parts) == 1 and relative.name.lower().startswith(('license', 'copyright')))
        selected |= relative == Path('ispc.html')
        if not selected:
            continue
        assert member.isfile() and member.size <= 1024 * 1024 * 512
        path = destination / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        assert not path.exists()
        data = package.extractfile(member).read()
        assert len(data) == member.size
        path.write_bytes(data)
        os.chmod(path, member.mode & 0o777)
        members[str(path.relative_to(root))] = {'archive_member': member.name,
            'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(), 'mode': member.mode & 0o777}
assert (destination / 'bin/ispc').is_file()
assert any('license' in Path(n).name.lower() for n in members)
archive.unlink()
record = {'scope': 'Pinned official compiler acquisition only. No compiler execution, source generation, shader/native runtime or acceptance.',
    'release_api': api, 'release_tag': tag, 'release_id': release['id'],
    'published_at': release['published_at'], 'asset': {k: asset[k] for k in ('id', 'name', 'size', 'digest', 'browser_download_url')},
    'extracted_members': members, 'verified_download_removed_after_extraction': True,
    'unused_GPU_runtime_and_examples_not_extracted': True}
(work / 'compiler-acquisition.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'compiler_acquired': True, 'download_bytes': total,
                  'retained_members': len(members), 'compiled_or_executed': False}))
