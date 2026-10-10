"""Bind complete generated payloads to the original exported PE consumers."""
from pathlib import Path
import hashlib, json, re, struct, subprocess, sys

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
BASELINE = 'eef7baf23f2e8c9906a467a6eef9ba7b6317d5ee'
CANDIDATE = 'bd82e523227a1b7223f43984af10291711925e44'
CHANGED = {'kTransportShader', 'kTransportFp64Shader',
           'kTransportPortableShader', 'kTransportPortableFp64Shader',
           'kTransportPortableNormalSumShader', 'kTransportFmaShader'}

def sha(data):
    return hashlib.sha256(data).hexdigest()

def arrays(data):
    result = {}
    pattern = r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};'
    for count, name, body in re.findall(pattern, data.decode(), re.S):
        words = [int(word[:-1], 0) for word in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u', body)]
        assert len(words) == int(count) and name not in result
        result[name] = struct.pack('<' + 'I' * len(words), *words)
    assert len(result) == 40
    return result

def sections(data):
    assert data[:2] == b'MZ'
    pe = struct.unpack_from('<I', data, 60)[0]
    assert data[pe:pe + 4] == b'PE\0\0'
    count = struct.unpack_from('<H', data, pe + 6)[0]
    optional = struct.unpack_from('<H', data, pe + 20)[0]
    result = []
    for i in range(count):
        offset = pe + 24 + optional + 40 * i
        virtual_size, virtual_address, size, start = struct.unpack_from('<IIII', data, offset + 8)
        flags = struct.unpack_from('<I', data, offset + 36)[0]
        result.append({'name': data[offset:offset + 8].rstrip(b'\0').decode('ascii'),
                       'raw_start': start, 'raw_bytes': size, 'virtual_bytes': virtual_size,
                       'virtual_address': virtual_address, 'characteristics': flags})
    return result

assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == CANDIDATE
assert not subprocess.check_output(['git', 'status', '--porcelain'])
label, = sys.argv[1:]
assert label == 'candidate'
revision = BASELINE if label == 'baseline' else CANDIDATE
header = WORK / (label + '-retained_kernels.h')
expected_header = ('3f7aca3372bb872934e978fb0270f2375ce8d2361555b35c8b891e42e07161cc'
                   if label == 'baseline' else
                   '8f5792677e201df1c7e409ca15afeaca07600fff570ec583b7ac41998380e963')
header_bytes = header.read_bytes()
assert sha(header_bytes) == expected_header
payloads = arrays(header_bytes)
baseline = arrays((WORK / 'baseline-retained_kernels.h').read_bytes())
changed = {name for name in payloads if payloads[name] != baseline[name]}
assert changed == (set() if label == 'baseline' else CHANGED)
accepted = json.loads((ROOT / 'attestations/native-vulkan/3045067/first-use-and-modern-original-controls/native-array-bindings.json').read_text())
for name, data in baseline.items():
    assert {'bytes': len(data), 'sha256': sha(data)} == {k: accepted['arrays'][name][k] for k in ('bytes', 'sha256')}
stage = ROOT / 'bin/windows-msvc' / ('native-' + revision[:7])
products, bindings = {}, {}
for name, relative in {'sirius.exe': 'src/sirius/app/Release/sirius.exe',
                       'sirius_backend_tests.exe': 'tests/backend/Release/sirius_backend_tests.exe'}.items():
    path = stage / relative
    data = path.read_bytes()
    sec = sections(data)
    products[name] = {'path': str(path.relative_to(ROOT)), 'bytes': len(data),
                      'sha256': sha(data), 'sections': sec}
    for array, payload in payloads.items():
        offsets, start = [], 0
        while (offset := data.find(payload, start)) != -1:
            offsets.append(offset)
            start = offset + 1
        assert offsets, (name, array)
        for offset in offsets:
            assert any(s['name'] == '.rdata' and s['characteristics'] & 0x40000000
                       and not s['characteristics'] & (0x80000000 | 0x20000000)
                       and s['raw_start'] <= offset
                       and offset + len(payload) <= s['raw_start'] + s['raw_bytes']
                       for s in sec), (name, array, offset)
        item = bindings.setdefault(array, {'bytes': len(payload), 'sha256': sha(payload),
                                          'consumer_bindings': {}})
        item['consumer_bindings'][name] = {'whole_payload_file_offsets': offsets, 'section': '.rdata'}
record = {'source_revision': revision, 'live_source_revision': CANDIDATE, 'label': label,
          'producer_header': str(header.relative_to(ROOT)), 'producer_header_sha256': sha(header_bytes),
          'changed_arrays': sorted(changed), 'unchanged_arrays': 40 - len(changed),
          'arrays': bindings, 'products': products, 'pass': True,
          'scope': 'Complete exact generated payloads in original official app/backend PE readable, nonwritable, nonexecutable .rdata. No all-array runtime selection, cross-provider identity, speed or scientific/release qualification claim.'}
output = WORK / (label + '-native-array-bindings.json')
assert not output.exists()
output.write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'label': label, 'arrays': len(bindings), 'changed': len(changed), 'products': len(products)}))
