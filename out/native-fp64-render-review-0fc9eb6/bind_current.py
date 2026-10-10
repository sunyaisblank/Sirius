"""Join current official products, whole arrays, inputs and original controls."""
from pathlib import Path
import hashlib
import json
import re
import struct
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
HEAD = json.loads((WORK / 'source.json').read_text())['source_revision']
STAGE = ROOT / 'bin/windows-msvc' / ('native-' + HEAD[:7])

def sha(b):
    return hashlib.sha256(b).hexdigest()

def identity(p):
    b = p.read_bytes()
    return {'bytes': len(b), 'sha256': sha(b)}

def write(name, value):
    (WORK / name).write_text(json.dumps(value, indent=2) + '\n')

assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == HEAD
assert not subprocess.check_output(['git', 'status', '--porcelain'])
header = ROOT / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
assert identity(header)['sha256'] == '3f7aca3372bb872934e978fb0270f2375ce8d2361555b35c8b891e42e07161cc'
arrays = {}
for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', header.read_text(), re.S):
    words = [int(w[:-1], 0) for w in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u', body)]
    assert len(words) == int(count) and name not in arrays
    arrays[name] = struct.pack('<' + 'I' * len(words), *words)
assert len(arrays) == 40
restored = json.loads((ROOT / 'attestations/native-vulkan/0fc9eb6/paired-transport-restoration/diagnostic/build-and-arrays.json').read_text())
for name, payload in arrays.items():
    assert {'bytes': len(payload), 'sha256': sha(payload)} == restored['arrays'][name]
gate = json.loads((STAGE / 'generated/sirius/native_build_gate.json').read_text())
assert gate['source'] == {'revision': HEAD, 'clean': True}
products = {}
joins = {}
for name in ('sirius', 'sirius_backend_tests', 'sirius_render_tests'):
    record = gate['tested_artifacts'][name]
    assert record['root'] == 'build'
    path = STAGE / record['path']
    data = path.read_bytes()
    assert identity(path) == {k: record[k] for k in ('bytes', 'sha256')}
    assert data[:2] == b'MZ'
    pe = struct.unpack_from('<I', data, 60)[0]
    assert data[pe:pe + 4] == b'PE\0\0'
    count = struct.unpack_from('<H', data, pe + 6)[0]
    optional = struct.unpack_from('<H', data, pe + 20)[0]
    sections = []
    for i in range(count):
        offset = pe + 24 + optional + 40 * i
        size, start = struct.unpack_from('<II', data, offset + 16)
        flags = struct.unpack_from('<I', data, offset + 36)[0]
        sections.append({'name': data[offset:offset + 8].rstrip(b'\0').decode(), 'start': start, 'bytes': size, 'flags': flags})
    products[name] = {'path': str(path.relative_to(ROOT)), **identity(path)}
    for array, payload in arrays.items():
        offsets, start = [], 0
        while (offset := data.find(payload, start)) != -1:
            offsets.append(offset)
            start = offset + 1
        assert offsets, (name, array)
        for offset in offsets:
            assert any(s['name'] == '.rdata' and s['flags'] & 0x40000000 and not s['flags'] & (0x80000000 | 0x20000000)
                       and s['start'] <= offset and offset + len(payload) <= s['start'] + s['bytes'] for s in sections)
        item = joins.setdefault(array, {'bytes': len(payload), 'sha256': sha(payload), 'consumers': {}})
        item['consumers'][name] = offsets
write('current-array-bindings.json', {'source_revision': HEAD, 'pass': True, 'header': identity(header), 'arrays': joins,
                                     'products': products, 'whole_payload_joins': 120, 'full_qualification_claimed': False})
files = {}
def bind(path):
    path = path.resolve()
    assert path.is_file() and not path.is_symlink()
    key = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
    record = identity(path)
    assert key not in files or files[key] == record
    files[key] = record

for rel in subprocess.check_output(['git', 'ls-files'], text=True).splitlines():
    p = ROOT / rel
    assert p.read_bytes() == subprocess.check_output(['git', 'show', HEAD + ':' + rel])
    bind(p)
for p in STAGE.rglob('*'):
    if p.is_file():
        bind(p)
for platform in ('windows', 'macos'):
    bundle = ROOT / 'attestations/native-build' / HEAD[:7] / (platform + '-build')
    for p in bundle.rglob('*'):
        if p.is_file():
            bind(p)
for name in ('source.json', 'plan.json', 'original-registration.json', 'historical-current-dopri-reuse.json', 'independent-prerequisite-review.json', 'prerequisite-comment.md',
             'prerequisite-comment-readback.json', 'current-array-bindings.json', 'export-ci-terminal.json',
             'export-jobs-terminal.json', 'export-artifacts.json', 'windows-export-download.json', 'macos-export-download.json',
             'reconstruction-windows-' + HEAD[:7] + '.json', 'download_exports.py', 'reconstruct_native.py',
             'bind_current.py', 'run_native_control.ps1', 'emergency_native_cleanup.ps1', 'run_current.py',
             'terminal-audit.ps1', 'run_terminal_current.py',
             'independent-owner-review.json', 'native-preregister-comment.md', 'native-preregister-readback.json'):
    bind(WORK / name)
bind(header)
providers = {
    '/mnt/c/Windows/System32/vulkan-1.dll': {'bytes': 1719408, 'sha256': '45ee2efc5c6d3986da8181775610f308122c752e55e703340d9733bbbc4b9f1a'},
    '/mnt/c/Windows/System32/DriverStore/FileRepository/amdvlk.inf_amd64_914ba89eaaafdf60/amdvlk64.dll': {'bytes': 103213064, 'sha256': '5c0561968155f4b0aa80fe9a5bbd17072d5942ecb906c687e749dc92058a8650'},
    '/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe': {'bytes': 495616, 'sha256': '8bb6fa8c283b4d92120b1ef249a9b311b0f804d4cabbe9981159976c8be76a5e'}
}
for path, record in providers.items():
    assert identity(Path(path)) == record
    bind(Path(path))
write('execution-bindings.json', {'source_revision': HEAD, 'stage': str(STAGE.relative_to(ROOT)),
                                  'gate_sha256': identity(STAGE / 'generated/sirius/native_build_gate.json')['sha256'],
                                  'files': files, 'providers': providers, 'whole_array_joins': 120,
                                  'full_qualification_claimed': False})
print(json.dumps({'source_revision': HEAD, 'bound_files': len(files), 'arrays': 40, 'whole_array_joins': 120, 'dispatches': 0}))
