from pathlib import Path
import hashlib
import json
import os
import shutil
import stat
import subprocess
import tempfile

root = Path.cwd()
parent = 'a75d4b547e4456b02b802ee2cee0f9e53974ef9d'
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == 'f5d649bf5067d1d68dc9215462f37830122f0540'
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
work = root / 'out/macos-moltenvk-identity-fixed'
payload = work / 'payload'
bundle = root / 'attestations/native-vulkan/4bcd987/macos-provider-identity-fixed'
assert not bundle.exists()
def encoded(x):
    return (json.dumps(x, indent=2, sort_keys=True) + '\n').encode()
def pin(p):
    raw = p.read_bytes()
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}
freeze_path = payload / 'input-freeze.json'
assert pin(freeze_path) == {'bytes': 2540081, 'sha256': 'a0b5986daccff1453e11361feb8b8a5a020a6e83c1943deec4da943f3e4ca8d4'}
freeze = json.loads(freeze_path.read_text())
actual_files, actual_links = {}, {}
for base, ds, fs in os.walk(payload, followlinks=False):
    for name in [*ds, *fs]:
        p = Path(base) / name
        relative = str(p.relative_to(payload))
        if relative == 'input-freeze.json':
            continue
        if p.is_symlink():
            assert p.resolve().is_relative_to(payload)
            actual_links[relative] = os.readlink(p)
        elif p.is_file():
            actual_files[relative] = {**pin(p), 'mode': stat.S_IMODE(p.stat().st_mode)}
assert actual_files == freeze['files'] and actual_links == freeze['symlinks']
assert len(actual_files) == 11803 and len(actual_links) == 10
content = json.loads((payload / 'content-id.json').read_text())
assert hashlib.sha256(encoded(content['preimage'])).hexdigest() == content['content_id'] == '76cb8b212bc049261146b526d75ebddafe017b18fe1f4c2c7065c3024c4ffe9e'
assert content['preimage']['build_profile'] == json.loads((payload / 'build-profile.json').read_text())
original = payload / 'historical/original'
assert pin(original / 'input-freeze.json') == {'bytes': 2531127, 'sha256': '71ca1547b3fb616ac99b9e8871a7e2410f25443131a6652b8b1f1c3464f526c9'}
reuse = json.loads((payload / 'translation-reuse.json').read_text())
old = json.loads((original / 'input-freeze.json').read_text())
changed = set(reuse['metadata_only_input_differences'])
assert changed == {'build-profile.json', 'content-id.json', 'patches/MoltenVK-diagnostic-identity.patch', 'sources/moltenvk-patched/MoltenVK/MoltenVK/GPUObjects/MVKDevice.mm'}
for name, row in old['files'].items():
    before = {k: row[k] for k in ('bytes', 'sha256')}
    p = payload / reuse['relocated_legacy_files'].get(name, name)
    if name in changed:
        assert pin(original / name) == before
        assert reuse['metadata_only_input_differences'][name] == {'original': before, 'current': pin(p)}
    else:
        assert pin(p) == before, name
device = 'sources/moltenvk-patched/MoltenVK/MoltenVK/GPUObjects/MVKDevice.mm'
old_statement = '\tstrcpy(supportedProps12.driverInfo, MVK_VERSION_STRING "+sirius-issue48-unforced." MVK_DIAGNOSTIC_CONTENT_ID);'
new_statement = '\tstd::snprintf(supportedProps12.driverInfo, sizeof(supportedProps12.driverInfo),\n\t\t\t\t  "%s+sirius-issue48-unforced.%s", MVK_VERSION_STRING, MVK_DIAGNOSTIC_CONTENT_ID);'
assert (original / device).read_text().replace('#include <cmath>\n', '#include <cmath>\n#include <cstdio>\n').replace(old_statement, new_statement) == (payload / device).read_text()
control = root / 'out/macos-provider-metadata-control'
for relative in ('original-compile', 'corrected-compile', 'runtime-owner/corrected-run'):
    o = json.loads((control / relative / 'owner.json').read_text())
    assert o['passed'] and o['source_unchanged'] and o['whole_inputs_unchanged']
    assert not o['stop_reason'] and not o['cleanup_errors'] and not o['remaining_owned_processes']
    assert all(b['absent'] for b in o['observed_births'])
assert b'expected' in (control / 'original-compile/stderr.log').read_bytes()
assert (control / 'runtime-owner/corrected-run/stdout.log').read_text() == content['driver_info_expected'] + '\n'
bundle.mkdir(parents=True)
shutil.copytree(payload, bundle / 'payload', symlinks=True)
shutil.copytree(control, bundle / 'metadata-control', symlinks=True)
shutil.copy2(root / 'bin/linux-gcc/macos-provider-metadata-control/corrected', bundle / 'metadata-control/corrected-elf')
for name in ('prepare_identity_fixed.py', 'prepared-payload.json', 'provider-source-join.json'):
    shutil.copy2(work / name, bundle / name)
shutil.copy2(root / 'out/prepare-macos-metadata-control.py', bundle / 'metadata-control/prepare-control.py')
shutil.copy2(Path(__file__), bundle / Path(__file__).name)
result = {'schema': 'sirius-corrected-diagnostic-provider-payload-v1', 'product_source': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
    'new_freeze': pin(freeze_path), 'content_id': content['content_id'], 'payload_files': len(actual_files), 'payload_links': len(actual_links),
    'review': 'Complete freeze, historical input map, exact four-file metadata/profile whitelist, canonical content ID and two-hunk source correction verified.',
    'compile_control': 'Original exact statement fails against actual pinned headers; corrected exact statement compiles and emits exact expected identity. Linux GCC control only.',
    'translation': 'Original exact translation reused with historical freeze; no fresh translation claimed.',
    'parent_runtime_action': 'Parent control corrected-run action was never used; separate runtime owner seals completed compiler output.',
    'native_provider_built_or_executed': False, 'qualification_claimed': False}
(bundle / 'result.json').write_bytes(encoded(result))
entries = {}
for p in sorted(bundle.rglob('*')):
    if p.is_symlink():
        entries[str(p.relative_to(bundle))] = {'type': 'symlink', 'target': os.readlink(p)}
    elif p.is_file():
        entries[str(p.relative_to(bundle))] = {'type': 'file', **pin(p), 'mode': stat.S_IMODE(p.stat().st_mode)}
(bundle / 'manifest.json').write_bytes(encoded({'files': entries, 'qualification_claimed': False}))
fd, index = tempfile.mkstemp(dir=root / 'out', prefix='preserve-identity-')
os.close(fd)
os.unlink(index)
env = dict(os.environ, GIT_INDEX_FILE=index)
try:
    subprocess.run(['git', 'read-tree', '--empty'], env=env, check=True)
    rows, verify = [], []
    for p in sorted(bundle.rglob('*')):
        if p.is_symlink():
            mode, raw = '120000', os.readlink(p).encode()
            blob = subprocess.check_output(['git', 'hash-object', '-w', '--stdin'], input=raw).decode().strip()
        elif p.is_file():
            mode, raw = ('100755' if p.stat().st_mode & 0o111 else '100644'), p.read_bytes()
            blob = subprocess.check_output(['git', 'hash-object', '-w', '--', str(p)]).decode().strip()
        else:
            continue
        rows.append(f'{mode} {blob}\t{p.relative_to(root)}\0')
        verify.append((str(p.relative_to(root)), hashlib.sha256(raw).hexdigest()))
    subprocess.run(['git', 'update-index', '-z', '--index-info'], input=''.join(rows).encode(), env=env, check=True)
    tree = subprocess.check_output(['git', 'write-tree'], env=env, text=True).strip()
    commit = subprocess.check_output(['git', 'commit-tree', tree, '-p', parent], input=b'Correct diagnostic provider metadata and preserve actual-header controls (#48)\n', env=env).decode().strip()
    tag = 'evidence/issue-48-provider-identity-fixed'
    subprocess.run(['git', 'tag', tag, commit], check=True)
    for name, digest in verify:
        assert hashlib.sha256(subprocess.check_output(['git', 'show', commit + ':' + name])).hexdigest() == digest
    payload_path = str((bundle / 'payload').relative_to(root))
    payload_tree = subprocess.check_output(['git', 'rev-parse', commit + ':' + payload_path], text=True).strip()
    print(json.dumps({'commit': commit, 'tag': tag, 'payload_path': payload_path, 'payload_tree': payload_tree, 'files': len(verify), 'all_git_bytes_verified': True, 'manifest': pin(bundle / 'manifest.json')}, indent=2))
finally:
    if Path(index).exists():
        Path(index).unlink()
