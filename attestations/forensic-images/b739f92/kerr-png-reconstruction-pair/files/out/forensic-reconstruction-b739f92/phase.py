"""Four finite preparation actions; deliberately contains no render action."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import time

sys.dont_write_bytecode = True
work = Path(__file__).resolve().parent
root = work.parent.parent
phase, folder_arg = sys.argv[1:]
folder = Path(folder_arg).resolve(strict=True)
assert folder.parent == work and folder.name == phase
manifest = json.loads((work / 'preparation.json').read_text())
source = work / 'source'
build = root / 'bin/linux-gcc/forensic-reconstruction-b739f92'
if phase.endswith('-scaffold'):
    build = build / 'scaffold-retry'
operation = phase.removesuffix('-scaffold')
assert phase in manifest['preparation_commands']
for name in ('CC', 'CXX', 'CFLAGS', 'CXXFLAGS', 'CPPFLAGS', 'LDFLAGS', 'CMAKE_TOOLCHAIN_FILE', 'CMAKE_PROJECT_INCLUDE', 'CMAKE_PROJECT_INCLUDE_BEFORE', 'CMAKE_PREFIX_PATH', 'CPATH', 'CPLUS_INCLUDE_PATH', 'LIBRARY_PATH', 'LD_LIBRARY_PATH', 'LD_PRELOAD', 'LD_AUDIT', 'GIT_CONFIG', 'GIT_CONFIG_GLOBAL', 'GIT_CONFIG_SYSTEM', 'GIT_CONFIG_COUNT'):
    assert not os.environ.get(name), ('undeclared build override', name)
env = os.environ.copy()
env.update(manifest['child_environment']['set'])
for name in manifest['child_environment']['unset']:
    env.pop(name, None)
records = []

def identity(path):
    data = path.read_bytes()
    return {'path': str(path), 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

def run(argv, *, expected=0):
    index = len(records)
    start = time.monotonic()
    result = subprocess.run(argv, cwd=source, env=env, capture_output=True)
    (folder / f'{index:02d}.stdout').write_bytes(result.stdout)
    (folder / f'{index:02d}.stderr').write_bytes(result.stderr)
    item = {'argv': argv, 'cwd': str(source), 'exit': result.returncode, 'expected_exit': expected,
            'wall_seconds': time.monotonic() - start}
    records.append(item)
    (folder / 'commands.json').write_text(json.dumps(records, indent=2) + '\n')
    assert result.returncode == expected, item
    return json.loads(result.stdout) if phase == 'comparator-controls' else result

def require_prior(name):
    owner = json.loads((work / name / 'owner.json').read_text())
    assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['remaining_owned_processes'] == [] and not owner['cleanup_errors'] and owner['stop_reason'] is None
    assert owner['observed_births'] and all(item['absent'] for item in owner['observed_births'])
    saved_receipt = json.loads((work / name / 'receipt.json').read_text())
    assert saved_receipt['passed'] and saved_receipt['phase'] == name
    for saved in saved_receipt['products']:
        assert identity(Path(saved['path'])) == saved, ('prior output changed', saved['path'])

if operation == 'configure-renderers':
    assert not (build / 'renderers').exists(), 'fresh private build required'
    (build / 'renderers').mkdir(parents=True)
    for argv in manifest['preparation_commands'][phase]:
        run(argv)
    cache = build / 'renderers/CMakeCache.txt'
    text = cache.read_text()
    assert f'CMAKE_HOME_DIRECTORY:INTERNAL={source}' in text
    assert 'SIRIUS_BUILD_LEGACY:BOOL=ON' in text and 'SIRIUS_MANDATORY_TESTS:BOOL=OFF' in text
    assert 'GLFW_INSTALL:BOOL=OFF' in text
    assert 'CMAKE_DISABLE_FIND_PACKAGE_Vulkan:UNINITIALIZED=TRUE' in text or 'CMAKE_DISABLE_FIND_PACKAGE_Vulkan:BOOL=TRUE' in text
    products = [cache, build / 'renderers/build.ninja', build / 'renderers/compile_commands.json']
elif operation == 'build-renderers':
    require_prior('configure-renderers-scaffold' if phase.endswith('-scaffold') else 'configure-renderers')
    for argv in manifest['preparation_commands'][phase]:
        run(argv)
    products = [build / 'renderers/bin/Sirius', build / 'renderers/bin/sirius',
                build / 'renderers/compile_commands.json', build / 'renderers/build.ninja']
elif operation == 'build-comparator':
    assert not (build / 'comparator').exists(), 'fresh private comparator build required'
    (build / 'comparator').mkdir(parents=True)
    for argv in manifest['preparation_commands'][phase]:
        run(argv)
    products = [build / 'comparator/image_compare', build / 'comparator/miniz.o']
else:
    require_prior('build-comparator')
    receipt = json.loads((work / 'build-comparator/receipt.json').read_text())
    for saved in receipt['products']:
        assert identity(Path(saved['path'])) == saved, 'compiled comparator changed'
    binary = str(build / 'comparator/image_compare')
    fixtures = folder / 'fixtures'
    assert run([binary, 'fixtures', str(fixtures)])['fixtures_created']
    before = {p.name: identity(p) for p in sorted(fixtures.iterdir())}
    for fmt in ('png', 'exr', 'ppm'):
        reference = str(fixtures / f'reference.{fmt}')
        common = [binary, 'compare', fmt, '3', '2', reference]
        same = run(common + [reference])
        assert same['equal'] and same['encoded_files_equal'] and same['mismatched_samples'] == 0
        metadata = run(common + [str(fixtures / f'metadata.{fmt}')])
        assert metadata['equal'] and not metadata['encoded_files_equal'] and metadata['mismatched_samples'] == 0
        changed = run(common + [str(fixtures / f'perturbed.{fmt}')], expected=1)
        assert not changed['equal'] and not changed['encoded_files_equal']
        assert changed['mismatched_pixels'] == changed['mismatched_samples'] == 1
        assert same['compared_pixels'] == metadata['compared_pixels'] == changed['compared_pixels'] == 6
        assert same['compared_samples'] == 6 * (3 if fmt == 'ppm' else 4)
        assert 'error' in run(common + [str(fixtures / 'missing')], expected=2)
        assert 'error' in run([binary, 'compare', fmt, '4', '2', reference, reference], expected=2)
    assert 'error' in run([binary, 'compare', 'png', '3', '2', str(fixtures / 'reference.exr'), str(fixtures / 'reference.exr')], expected=2)
    truncated = folder / 'truncated.png'
    truncated.write_bytes((fixtures / 'reference.png').read_bytes()[:-1])
    assert 'error' in run([binary, 'compare', 'png', '3', '2', str(truncated), str(truncated)], expected=2)
    assert before == {p.name: identity(p) for p in sorted(fixtures.iterdir())}
    products = sorted(fixtures.iterdir()) + [truncated]

receipt = {'evidence_kind': 'NEW_RECONSTRUCTION', 'phase': phase, 'passed': True,
           'commands': records, 'child_environment': manifest['child_environment'],
           'products': [identity(p) for p in products], 'historical_tape_identity_claim': False}
(folder / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({'phase': phase, 'passed': True, 'products': len(products), 'rendered': False}), flush=True)
