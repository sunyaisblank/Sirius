"""Bind complete baseline40/candidate42 payloads in actual official PE consumers."""
from pathlib import Path
import hashlib, importlib.util, json, re, struct, subprocess, sys
sys.dont_write_bytecode = True
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
label, = sys.argv[1:]
assert label in ('baseline', 'candidate')
source = json.loads((WORK / 'source.json').read_text())
revision = source['baseline_revision'] if label == 'baseline' else source['source_revision']
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == source['source_revision']
assert not subprocess.check_output(['git', 'status', '--porcelain'])
header = ROOT / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
build = json.loads((WORK / 'build-and-arrays.json').read_text())
assert hashlib.sha256(header.read_bytes()).hexdigest() == build['whole_header_sha256']
baseline = json.loads((WORK / 'baseline-arrays.json').read_text())['arrays']
payloads = {}
for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', header.read_text(), re.S):
    words = [int(q[:-1], 0) for q in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u', body)]
    assert len(words) == int(count) and name not in payloads
    payloads[name] = struct.pack('<' + 'I' * len(words), *words)
assert len(payloads) == 42 and len(baseline) == 40
assert set(payloads) - set(baseline) == {'kDenseFmaShader', 'kDopriPhaseFmaShader'}
for name, payload in payloads.items():
    identity = {'bytes': len(payload), 'sha256': hashlib.sha256(payload).hexdigest()}
    assert identity == build['arrays'][name]
    if name in baseline:
        assert identity == {k: baseline[name][k] for k in ('bytes', 'sha256')}
if label == 'baseline':
    payloads = {name: payloads[name] for name in baseline}
stage = ROOT / 'bin/windows-msvc' / ('native-' + revision[:7])
spec = importlib.util.spec_from_file_location('attestation', ROOT / 'scripts/verify-attestation.py')
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)
g = v.load_build_gate_verifier()
gate_path = stage / 'generated/sirius/native_build_gate.json'
gate = g.validate_native_build_document(json.loads(gate_path.read_text()))
assert gate['source'] == {'revision': revision, 'clean': True}
g.verify_recorded_files(gate, stage / 'recorded-source', stage)
assert len(v.copy_qualification_test_inputs(gate_path, stage)) == 20
products, joins = {}, {}
for name in ('sirius', 'sirius_backend_tests', 'sirius_render_tests'):
    record = gate['tested_artifacts'][name]
    assert record['root'] == 'build'
    path = stage / record['path']
    data = path.read_bytes()
    identity = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
    assert identity == {k: record[k] for k in ('bytes', 'sha256')}
    assert data[:2] == b'MZ'
    pe = struct.unpack_from('<I', data, 60)[0]
    assert data[pe:pe + 4] == b'PE\0\0'
    count, optional = struct.unpack_from('<H', data, pe + 6)[0], struct.unpack_from('<H', data, pe + 20)[0]
    sections = []
    for i in range(count):
        at = pe + 24 + optional + 40 * i
        size, start = struct.unpack_from('<II', data, at + 16)
        flags = struct.unpack_from('<I', data, at + 36)[0]
        sections.append({'name': data[at:at + 8].rstrip(b'\0').decode(), 'start': start, 'bytes': size, 'flags': flags})
    products[name] = {'path': str(path.relative_to(ROOT)), **identity}
    for array, payload in payloads.items():
        offsets, start = [], 0
        while (offset := data.find(payload, start)) != -1:
            offsets.append(offset)
            start = offset + 1
        assert offsets, (label, name, array)
        for offset in offsets:
            assert any(s['name'] == '.rdata' and s['flags'] & 0x40000000 and not s['flags'] & (0x80000000 | 0x20000000)
                       and s['start'] <= offset and offset + len(payload) <= s['start'] + s['bytes'] for s in sections)
        item = joins.setdefault(array, {'bytes': len(payload), 'sha256': hashlib.sha256(payload).hexdigest(), 'consumers': {}})
        item['consumers'][name] = offsets
result = {'source_revision': revision, 'pass': True, 'label': label, 'arrays': joins, 'products': products,
          'original40_unchanged': True, 'whole_payload_joins': len(payloads) * 3,
          'scope': 'Complete actual readonly PE payload and20 input joins, no loaded-stage/execution/performance/fullqualification claim.'}
output = WORK / (label + '-native-array-bindings.json')
assert not output.exists()
output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'label': label, 'revision': revision, 'arrays': len(payloads), 'complete_readonly_joins': len(payloads) * 3}))
