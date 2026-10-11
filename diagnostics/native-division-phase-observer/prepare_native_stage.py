"""Validate the existing CI export and reconstruct its exact ignored native stage."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, importlib.util, json, re, shutil, struct, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
HEAD = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert HEAD == '451e2cbf833b87f51b1819421d68163a4b5fd1e5'
assert HEAD == subprocess.check_output(['git', 'rev-parse', '@{u}'], text=True).strip()
assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
ARTIFACT = WORK / 'artifact'
STAGE = ROOT / 'bin/windows-msvc/native-451e2cb-division-phase-observer'
assert not STAGE.exists() and not (WORK/'execution-plan.json').exists()
assert ROOT.resolve() == Path('/home/astra/.project/Sirius')
assert WORK == ROOT/'out/native-division-phase-observer'

def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))

required_reviews = ('host-source-review.json', 'reader-source-review.json', 'compiled-review.json', 'actual-placement-review.json')
review_paths = {Path(__file__), WORK/'native-preparation-review.json', WORK/'candidate-source.json', WORK/'baseline-source.json', WORK/'source-build.json', WORK/'ci-dispatch.json', WORK/'ci-completed.json'}
preparer_review = read(WORK/'native-preparation-review.json')
assert preparer_review['accepted'] and preparer_review['source_revision'] == HEAD
assert preparer_review['reviewed_seals'][str(Path(__file__))] == seal(Path(__file__))
for path, expected in preparer_review['reviewed_seals'].items():
    assert seal(Path(path)) == expected, path
    review_paths.add(Path(path))
for name in required_reviews:
    path = WORK/name
    review = read(path)
    assert review['accepted']
    if 'source_revision' in review: assert review['source_revision'] == HEAD
    assert str(path) in preparer_review['reviewed_seals'], name
    review_paths.add(path)
    for source, expected in review['reviewed_seals'].items():
        assert seal(Path(source)) == expected, source
        review_paths.add(Path(source))
build = read(WORK/'source-build.json')
assert build['accepted'] and build['source_revision'] == HEAD
for path, expected in build['input_seals'].items():
    assert seal(Path(path)) == expected, path
    review_paths.add(Path(path))
candidate = read(WORK/'candidate-source.json')
assert candidate['source_revision'] == HEAD and candidate['committed'] and candidate['pushed']
assert candidate['baseline_revision'] == '2053ef79c90752d37b0075b8137cfb9b571b6820'
for path, expected in candidate['changed_source_seals'].items(): assert seal(ROOT/path) == expected
assert subprocess.check_output(['git','rev-parse','HEAD^'],text=True).strip() == candidate['baseline_revision']
assert set(subprocess.check_output(['git','diff','--name-only','HEAD^','HEAD'],text=True).splitlines()) == set(candidate['changed_source_seals'])
ci = read(WORK/'ci-completed.json'); dispatch = read(WORK/'ci-dispatch.json')
assert ci['databaseId'] == dispatch['run_id'] == 38168463363 and ci['headSha'] == HEAD
assert ci['headBranch'] == candidate['development_branch'] and dispatch['inputs'] == {'full_qualification': False}
assert ci['status'] == 'completed' and ci['conclusion'] == 'success'
assert {j['name'] for j in ci['jobs'] if j['conclusion'] == 'success'} == {'integration-no-render','integration-macos-no-render','integration-windows-no-render'}
assert {j['name'] for j in ci['jobs'] if j['conclusion'] == 'skipped'} == {'linux-gate','linux-sanitizers','windows-build','macos-build'}
source_seals = {str(path):seal(path) for path in sorted(review_paths)}

def safe(root, relative):
    path = root / relative
    assert not Path(relative).is_absolute() and '..' not in Path(relative).parts
    assert path.resolve().is_relative_to(root.resolve())
    return path

document = json.loads((ARTIFACT / 'windows-build.json').read_text())
assert document['status'] == 'pass' and document['source_revision'] == HEAD
assert document['domains'] == ['windows-native-build']
assert document['platform'] == {'os': 'windows', 'wsl2': False}
assert document['claims']['compiled_topology_complete']
assert document['claims']['runtime_ready'] is False
payloads = {}
for name, record in document['artifacts'].items():
    path = safe(ARTIFACT, record['path'])
    assert seal(path) == {k: record[k] for k in ('bytes', 'sha256')}
    payloads.setdefault((record['bytes'], record['sha256']), path)

sys.path.insert(0, str(ROOT / 'scripts'))
spec = importlib.util.spec_from_file_location('native_gate_verifier', ROOT / 'scripts/verify-build-gate.py')
verifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)
gate = json.loads((ARTIFACT / 'native_build_gate.json').read_text())
verifier.validate_native_build_document(gate)
assert gate['source'] == {'revision': HEAD, 'clean': True}
STAGE.mkdir(parents=True)

def copy_record(record, destination):
    source = payloads[(record['bytes'], record['sha256'])]
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        assert destination.read_bytes() == source.read_bytes()
    else:
        shutil.copyfile(source, destination)
    assert seal(destination) == {k: record[k] for k in ('bytes', 'sha256')}

tested = {}
for name, record in gate['tested_artifacts'].items():
    assert record['root'] == 'build'
    tested[name] = safe(STAGE, record['path'])
    copy_record(record, tested[name])
for record in list(gate['product_artifacts'].values()) + list(gate['test_input_artifacts'].values()):
    if record['root'] == 'build':
        copy_record(record, safe(STAGE, record['path']))
    else:
        assert record['root'] == 'source'
        assert seal(safe(ROOT, record['path'])) == {k: record[k] for k in ('bytes', 'sha256')}

for name, record in gate['product_artifacts'].items():
    if name == 'sirius': continue
    relative = verifier.INSTALLED_PRODUCTS[name].removeprefix('share/sirius/')
    copy_record(record, safe(tested['sirius'].parent / 'resources', relative))
for consumer in ('sirius_backend_tests', 'sirius_render_tests'):
    resource = tested[consumer].parent / 'resources'
    for name, record in gate['test_input_artifacts'].items():
        copy_record(record, safe(resource, verifier.TEST_INPUT_PATHS[name]))
    for name in ('trace_spv', 'trace_fp32comp_spv', 'trace_fp64_spv'):
        copy_record(gate['product_artifacts'][name], safe(resource, verifier.INSTALLED_PRODUCTS[name].removeprefix('share/sirius/')))
verifier.verify_test_input_copies(tested, gate['product_artifacts'], gate['test_input_artifacts'])
gate_path = STAGE / 'generated/sirius/native_build_gate.json'
gate_path.parent.mkdir(parents=True, exist_ok=True)
shutil.copyfile(ARTIFACT / 'native_build_gate.json', gate_path)
exe = tested['sirius_backend_tests'].read_bytes()
# Bind complete payloads to bytes mapped as initialized read-only PE data.
assert exe[:2] == b'MZ'
pe = struct.unpack_from('<I',exe,0x3c)[0]
assert exe[pe:pe+4] == b'PE\0\0'
machine, count = struct.unpack_from('<HH',exe,pe+4)
optional_bytes = struct.unpack_from('<H',exe,pe+20)[0]
assert machine == 0x8664 and 0 < count <= 96 and optional_bytes >= 112
assert struct.unpack_from('<H',exe,pe+24)[0] == 0x20b
section_table = pe+24+optional_bytes
assert section_table+count*40 <= len(exe)
readonly = []
for i in range(count):
    base = section_table+i*40
    name = exe[base:base+8].rstrip(b'\0').decode('ascii')
    virtual_size, rva, raw_bytes, file_offset = struct.unpack_from('<IIII',exe,base+8)
    flags = struct.unpack_from('<I',exe,base+36)[0]
    assert not raw_bytes or file_offset+raw_bytes <= len(exe)
    if flags & 0x40 and flags & 0x40000000 and not flags & (0x80000000|0x20000000):
        readonly.append({'section':name,'rva':rva,'file_offset':file_offset,'mapped_bytes':min(virtual_size,raw_bytes),'characteristics':flags})
assert readonly

def payload_binding(data):
    matches = []
    for section in readonly:
        begin = section['file_offset']; body = exe[begin:begin+section['mapped_bytes']]
        position = body.find(data)
        while position >= 0:
            matches.append({**section,'payload_file_offset':begin+position,'payload_rva':section['rva']+position})
            position = body.find(data,position+1)
    assert matches, 'Exact full payload absent from initialized readonly PE sections'
    return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'readonly_PE_bindings':matches}

header = ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
assert seal(header) == read(WORK/'baseline-source.json')['production_header']
production = {}
for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',header.read_text(),re.S):
    words = [int(word[:-1],0) for word in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body)]
    assert len(words) == int(count) and name not in production
    production[name] = payload_binding(struct.pack('<'+'I'*len(words),*words))
assert len(production) == 40
embedded = {}
for name in ('retained_transport.spv', 'retained_endpoint.spv'):
    path = (WORK/'phase/retained_transport.spv' if name == 'retained_transport.spv' else ROOT/'attestations/native-vulkan/2053ef7/native-shader-clock-admission/files/phase/retained_endpoint.spv')
    expected = read(WORK/'host-source-review.json')['transport_observer_module' if name == 'retained_transport.spv' else 'endpoint_observer_module']
    assert seal(path) == expected
    embedded[name] = payload_binding(path.read_bytes())
    review_paths.add(path)
assert {path:seal(Path(path)) for path in sorted(source_seals)} == source_seals


providers = json.loads((ROOT / 'attestations/native-vulkan/2053ef7/native-shader-clock-placement/files/transport/plan.json').read_text())['providers']
for provider in providers:
    assert seal(Path(provider['linux_path'])) == {k: provider[k] for k in ('bytes', 'sha256')}
review_paths.update((ARTIFACT/'windows-build.json', ARTIFACT/'native_build_gate.json', ROOT/'scripts/verify-build-gate.py', ROOT/'attestations/native-vulkan/2053ef7/native-shader-clock-placement/files/transport/plan.json'))
inputs = {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p):seal(p) for p in sorted(set(p for p in STAGE.rglob('*') if p.is_file()) | review_paths)}
assert {path:seal(Path(path)) for path in sorted(source_seals)} == source_seals
assert HEAD == subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
assert HEAD == subprocess.check_output(['git','rev-parse','@{u}'],text=True).strip()
assert not subprocess.check_output(['git','status','--porcelain'],text=True)

plan = {
    'source_revision': HEAD, 'build_root': str(STAGE.relative_to(ROOT)),
    'executable_directory': str(tested['sirius_backend_tests'].parent.relative_to(ROOT)),
    'gate_sha256': seal(gate_path)['sha256'], 'input_seals': inputs,
    'providers': providers, 'exact_observer_modules': embedded,
    'whole_production_arrays_readonly_PE': production, 'production_header':seal(header),
    'exported_backend_executable':seal(tested['sirius_backend_tests']), 'ci_run_id':ci['databaseId'],
    'cases': {'device-timestamps': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults',
              'mixed-transport': 'RetainedComputeTest.MixedIndependentTransportLayersPreserveWordsAndRefusal',
              'factored-endpoints': 'RetainedComputeTest.FactoredEndpointsPreserveIndependentRootsAndBoundaryAdmission'},
    'purpose': 'Original numerical/ordered-stream/refusal/recovery controls with uncalibrated native-default shader-clock phase/division sidecars and exact input snapshots; original inputs/modes/assertions/8MiB/governors intact. Existing 320s/4GiB native tool guard, 400s bridge. No full science/release or ordinary-benefit claim.'}
(WORK / 'execution-plan.json').write_text(json.dumps(plan, indent=2) + '\n')
print(json.dumps({'source_revision': HEAD, 'stage_files': sum(p.is_file() for p in STAGE.rglob('*')), 'readonly_production_arrays':len(production), 'native_backend': seal(tested['sirius_backend_tests']), 'observer_modules': embedded}))
