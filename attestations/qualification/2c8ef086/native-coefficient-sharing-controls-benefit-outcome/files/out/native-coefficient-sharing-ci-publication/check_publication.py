import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
assert len(sys.argv) == 5, 'STAGE ACTUAL_PRIVATE_SHA ACTUAL_RUN_ID COMPLETED_RUN_JSON'
STAGE = ROOT / sys.argv[1]
REVISION = sys.argv[2]
RUN = int(sys.argv[3])
RUN_PATH = ROOT / sys.argv[4]
BASE = '2c8ef08693d46f57d5dea14efdbddd4874dac684'
assert STAGE.is_relative_to(ROOT / 'bin/windows-msvc') and STAGE.resolve(strict=True) == STAGE
assert re.fullmatch('[0-9a-f]{40}', REVISION) and RUN > 0
assert not (WORK / 'root-publication-check.json').exists()
for arm in ['baseline', 'candidate']:
    for suffix in ['stdout', 'stderr']:
        assert not (WORK / (arm + '-PE-imports.' + suffix)).exists()


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def norm(path):
    return str(path).replace('\\', '/').lower()


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


producer = read(STAGE / 'producer-result.json')
assert producer['schema'] == 'sirius-native-packed-matched-ci-producer-v1'
assert producer['producer_profile'] == 'native-coefficient-sharing-one-shader-only'
assert producer['canonical_CXX_source_revision'] == BASE
assert type(producer['private_compile_object_count']) is int and producer['private_compile_object_count'] == 5
assert type(producer['expected_private_compile_objects']) is int and producer['expected_private_compile_objects'] == 5
assert producer['source_revision'] == producer['github_sha'] == REVISION
assert producer['run_id'] == str(RUN) and str(producer['run_attempt']) == '1'
before = read(STAGE / 'build-inputs-before.json')
assert before == read(STAGE / 'build-inputs-after.json')
source = read(STAGE / 'source-before-configure.json')
prefix = norm(producer['source_root']) + '/'
source_rows = {}
for name, row in source.items():
    normalized = norm(name)
    assert normalized.startswith(prefix)
    relative = normalized[len(prefix):]
    assert relative not in source_rows
    source_rows[relative] = row
tree = subprocess.check_output(['git', 'ls-tree', '-rz', REVISION], cwd=ROOT)
expected = {}
objects = {}
for item in tree.split(b'\0'):
    if not item:
        continue
    metadata, name = item.split(b'\t', 1)
    mode, kind, blob = metadata.decode().split()
    assert kind == 'blob'
    if blob not in objects:
        data = subprocess.check_output(['git', 'cat-file', 'blob', blob], cwd=ROOT)
        objects[blob] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
    expected[norm(name.decode())] = objects[blob]
# The private source commit must conserve every canonical2c8 blob.
canonical_tree = subprocess.check_output(['git', 'ls-tree', '-rz', BASE], cwd=ROOT)
canonical_count = 0
for item in canonical_tree.split(b'\0'):
    if not item: continue
    metadata, name = item.split(b'\t', 1)
    mode, kind, oid = metadata.decode().split()
    assert kind == 'blob'
    if oid not in objects:
        data = subprocess.check_output(['git', 'cat-file', 'blob', oid], cwd=ROOT)
        objects[oid] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
    assert expected[norm(name.decode())] == objects[oid], name
    canonical_count += 1
assert set(expected) == set(source_rows)
for name, identity in expected.items():
    assert all(source_rows[name][key] == identity[key] for key in ['bytes', 'sha256']), name
build_rows = {norm(name): row for name, row in before.items()}
assert len(build_rows) == len(before)
tools = read(STAGE / 'selected-tools.json')
tool_rows = {}
for name, path in tools.items():
    row = build_rows[norm(path)]
    assert row['readable_stream_bytes'] > 0
    if row['bytes'] == 0:
        assert row['configured_alias_metadata'] is True and name == 'python'
        sibling = build_rows[norm(tools['python_nonempty_sibling'])]
        assert sibling['bytes'] > 0 and sibling['sha256'] == row['sha256']
    tool_rows[name] = row
run = read(RUN_PATH)
assert run['head_sha'] == REVISION and run['status'] == 'completed' and run['conclusion'] == 'success'
assert run['event'] == 'push' and run['run_attempt'] == 1 and run['id'] == RUN
assert run['path'] == '.github/workflows/native-coefficient-sharing-diagnostic.yml'
packaging = read(STAGE / 'runtime-packaging.json')
imports = {}
for arm in ['baseline', 'candidate']:
    folder = STAGE / 'payload' / arm
    payload_names = {item.name.lower(): item for item in folder.iterdir() if item.is_file()}
    assert len(payload_names) == len([item for item in folder.iterdir() if item.is_file()])
    for name, origin in packaging['runtime_origins'].items():
        assert pin(folder / name) == {key: origin[key] for key in ['bytes', 'sha256']}
    exe = folder / ('original_' + arm + '_controls.exe')
    result = subprocess.run(['/usr/bin/objdump', '-p', str(exe)], cwd=ROOT,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    (WORK / (arm + '-PE-imports.stdout')).write_bytes(result.stdout)
    (WORK / (arm + '-PE-imports.stderr')).write_bytes(result.stderr)
    names = re.findall(r'DLL Name:\s*(\S+)', result.stdout.decode())
    assert names and len(names) == len(set(names))
    assert 'vulkan-1.dll' in {name.lower() for name in names}
    assert not any(name.lower().endswith('140d.dll') or name.lower() == 'ucrtbased.dll' for name in names)
    for name in names:
        if name.lower().startswith(('vulkan-', 'msvcp', 'vcruntime', 'concrt')):
            assert name.lower() in payload_names, name
    imports[arm] = names
assert imports['baseline'] == imports['candidate']
log = (STAGE / 'build.stdout.log').read_text(encoding='utf-8-sig')
warnings = re.findall(r'\bwarning\s+([A-Z]\d+)\b', log, re.I)
assert not any(code.upper().startswith('LNK') for code in warnings)
assert len(tools) == 15
record = {'schema': 'sirius-native-coefficient-sharing-actual-publication-check-v1', 'passed': True,
          'producer_source_revision': REVISION, 'all_Git_source_files': len(expected),
          'all_source_bytes_match_immutable_Git': True, 'canonical2c8_files_conserved': canonical_count, 'before_after_build_rows': len(before),
          'selected_tools': tool_rows, 'PE_imports': imports,
          'runtime_packaging_versions': [packaging['msvc_version'], packaging['release_crt_version']],
          'warning_codes_including_summaries': {code: warnings.count(code) for code in sorted(set(warnings))},
          'actual_completed_one_attempt_CI_exact': True,
          'runtime_admitted': False,
          'scope': 'Actual source/API/build-input/selected-tool/staged SDK+CRT/import joins. Tool+adjacent-DLL/configured/fetched closure, not complete host/SDK headers/transitive DSOs. Static private object/PE/MAP and concrete phase review remain separate.'}
(WORK / 'root-publication-check.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({key: value for key, value in record.items() if key not in ['selected_tools', 'PE_imports']}))
