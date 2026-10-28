"""Acquire the published matching LLVM SDK inside the existing preset cache."""
from pathlib import Path, PurePosixPath
import copy
import hashlib
import json
import re
import subprocess
import sys
import tarfile

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-strict-ir-linking'
D = B / 'toolchain'
ASSET = 'LLVM-23.1.2-Linux-X64.tar.xz'
URL = 'https://github.com/llvm/llvm-project/releases/download/llvmorg-23.1.2/' + ASSET
SHA = 'b5ed9675149cc837c282e9b6962c276c9fa62863d5b2f91537b60848552995b7'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


phase, = sys.argv[1:]
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    owner = json.loads((W / 'emit/owner.json').read_text())
    assert owner['status'] == 'terminal' and owner['passed']
    assert not D.exists() and not (B / ASSET).exists() and not (B / (ASSET + '.part')).exists()
    actions = json.loads((W / 'actions.json').read_text())
    assert set(actions) == {'emit'}
    source = dict(version='23.1.2', url=URL, sha256=SHA,
        publisher='https://github.com/llvm/llvm-project/releases/expanded_assets/llvmorg-23.1.2',
        format_choice='Published xz asset uses the already installed Python lzma decoder; no new decompression dependency.',
        signature_status='Publisher signature/attestation are available; this acquisition verifies the published SHA256, not a signature.')
    write(W / 'toolchain-source.json', source)
    inputs = [*W.glob('*.py'), W / 'toolchain-source.json',
              R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
              R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py',
              Path('/usr/bin/curl').resolve(), Path('/usr/bin/ldd').resolve()]
    actions['acquire'] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=['/usr/bin/python3', '-B', str(W / 'acquire.py'), 'acquire'],
        frozen_inputs=[str(p.relative_to(R)) if p.is_relative_to(R) else str(p) for p in inputs],
        scope='One pinned official LLVM23.1.2 SDK acquisition, safe extraction and version/dependency check inside bin/linux-gcc. All setup cost included; no compiler study, native workload, math/FP flag change, global installation or runtime qualification.')
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared='acquire', executed=False, checksum=SHA)))
elif phase == 'acquire':
    archive, partial = B / ASSET, B / (ASSET + '.part')
    assert not D.exists() and not archive.exists()
    with partial.open('xb') as stream:
        subprocess.run(['/usr/bin/curl', '--fail', '--location', '--silent', '--show-error',
                        '--retry', '0', '--max-time', '240', '--proto', '=https', URL],
                       stdout=stream, check=True)
    assert pin(partial)['sha256'] == SHA, 'Published checksum mismatch; do not extract'
    partial.rename(archive)
    D.mkdir()
    root = None
    count = total = 0
    with tarfile.open(archive, mode='r|xz') as bundle:
        for original in bundle:
            parts = PurePosixPath(original.name).parts
            assert parts and not PurePosixPath(original.name).is_absolute() and '..' not in parts
            root = root or parts[0]
            assert parts[0] == root and (original.isfile() or original.isdir() or original.issym() or original.islnk())
            if len(parts) == 1:
                assert original.isdir()
                continue
            member = copy.copy(original)
            member.name = str(PurePosixPath(*parts[1:]))
            if member.islnk():
                link = PurePosixPath(member.linkname)
                assert not link.is_absolute() and '..' not in link.parts and link.parts[0] == root
                member.linkname = str(PurePosixPath(*link.parts[1:]))
            bundle.extract(member, path=D, filter='data')
            count += 1
            total += original.size if original.isfile() else 0
    tools, libraries = {}, {}
    for name in ('clang++', 'llvm-link', 'opt', 'llc', 'llvm-as', 'llvm-dis', 'llvm-config'):
        path = D / 'bin' / name
        real = path.resolve(strict=True)
        assert real.is_relative_to(D) and real.is_file()
        version = subprocess.check_output([str(path), '--version'], text=True, stderr=subprocess.STDOUT)
        assert '23.1.2' in version, (name, version)
        resolution = subprocess.check_output(['/usr/bin/ldd', '-r', str(real)], text=True, stderr=subprocess.STDOUT)
        assert 'not found' not in resolution and 'undefined symbol' not in resolution, (name, resolution)
        for value in re.findall(r'(/[^\s()]+)', resolution):
            library = Path(value)
            if library.is_file():
                library = library.resolve(strict=True)
                libraries[str(library)] = pin(library)
        tools[name] = dict(path=str(path.relative_to(R)), resolved=str(real.relative_to(R)),
                           **pin(real), version=version, dependency_resolution=resolution)
    prefix = subprocess.check_output([str(D / 'bin/llvm-config'), '--prefix'], text=True).strip()
    assert Path(prefix).resolve(strict=True) == D.resolve(strict=True)
    write(W / 'toolchain.json', dict(acquired=True, archive=dict(path=str(archive.relative_to(R)), **pin(archive)),
        extracted_root=root, extracted_members=count, unpacked_regular_bytes=total,
        tools=tools, runtime_libraries=libraries, relocated_prefix_correct=True,
        scope='Published SHA256, exact tool version, contained extraction and actual dynamic dependencies verified. No signature claim, source compilation, constrained conversion/inlining, numerical or speed/adoption qualification. SDK cache remains local to this checkout.'))
    print(json.dumps(dict(acquired=True, version='23.1.2', verified_tools=list(tools), members=count, unpacked_bytes=total)))
else:
    raise ValueError(phase)
