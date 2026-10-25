"""One copied-object rename/fixture build or one full-word raw ABI comparison."""
from pathlib import Path
import hashlib
import json
import os
import re
import shlex
import struct
import subprocess
import sys
from elf_proof import pin, rename_proof, associate

R = Path.cwd()
W = Path(__file__).resolve().parent


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def unchanged(prepared):
    for name, expected in prepared['source_inputs'].items():
        assert pin(R / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name


def identity(path):
    return {'path': str(path.relative_to(R)), **pin(path)}


mode, output = sys.argv[1:]
folder = Path(output); assert mode in ('compile', 'run') and folder == W / mode
assert not any(os.environ.get(k) for k in ('CPATH', 'CPLUS_INCLUDE_PATH', 'LIBRARY_PATH', 'GCC_EXEC_PREFIX', 'COMPILER_PATH'))
prepared = json.loads((W / 'prepared.json').read_text()); unchanged(prepared)
baseline, driver, executable = [R / prepared[k] for k in ('baseline', 'driver', 'executable')]
if mode == 'compile':
    result_path = W / 'compile-result.json'; assert not result_path.exists()
    record = {'passed': False, 'prepared': pin(W / 'prepared.json'), 'commands': [],
              'all8_native_objects': prepared['objects'], 'candidate_transport_object': prepared['candidate_transport_object'],
              'original_baseline_object': prepared['original_baseline_object'],
              'new_baseline_math_compilations': 0, 'diagnostic_CPP_compilations': 0, 'diagnostic_links': 0}
    write(result_path, record)
    old = R / prepared['original_baseline_object']['path']
    globals_text = subprocess.check_output(['/usr/bin/nm', '-g', '--defined-only', str(old)], text=True)
    (folder / 'original-global-definitions.txt').write_text(globals_text)
    assert {line.split()[-1] for line in globals_text.splitlines() if line.split()} == set(prepared['rename_mapping'])
    assert not baseline.exists() and not driver.exists() and not executable.exists()
    commands = [prepared['rename_argv'], prepared['compile_argv'], prepared['link_argv']]
    for index, argv in enumerate(commands):
        code = subprocess.run(argv, cwd=R).returncode
        record['commands'].append({'argv': argv, 'cwd': str(R), 'exit': code}); write(result_path, record)
        if code:
            raise SystemExit(code)
        if index == 0:
            proof = rename_proof(old, baseline, prepared['rename_mapping'])
            assert pin(old) == {k: prepared['original_baseline_object'][k] for k in ('bytes', 'sha256')}
            write(W / 'rename-proof.json', proof); record['renamed_baseline_object'] = identity(baseline)
        if index == 1:
            record['diagnostic_CPP_compilations'] = 1; record['driver_object'] = identity(driver)
            text = (folder / 'driver.d').read_text().replace('\\\n', ' ')
            dependencies = {}
            for name in shlex.split(text.split(':', 1)[1]):
                p = Path(name)
                p = (p if p.is_absolute() else R / p).resolve(strict=True)
                key = str(p.relative_to(R)) if p.is_relative_to(R) else str(p)
                expected = prepared['source_inputs'].get(key, prepared['source_inputs'].get(str(p)))
                assert expected is not None and pin(p) == {k: expected[k] for k in ('bytes', 'sha256')}, key
                dependencies[key] = pin(p)
            record['actual_compiled_dependencies'] = dependencies
        if index == 2:
            record['diagnostic_links'] = 1
        write(result_path, record)
    new = prepared['candidate_transport_object']
    new_absolute = {**new, 'path': str(R / new['path'])}
    old_absolute = {'path': str(baseline), **pin(baseline)}
    proof = associate(executable, [old_absolute, new_absolute])
    # Standard checkout-relative receipt names; association itself used exact absolute files.
    proof['executable'] = identity(executable)
    proof['objects'][0]['object'] = identity(baseline)
    proof['objects'][1]['object'] = new
    write(W / 'link-proof.json', proof)
    libraries_text = subprocess.check_output(['/usr/bin/ldd', '-r', str(executable)], text=True, stderr=subprocess.STDOUT)
    (folder / 'runtime-resolution.txt').write_text(libraries_text)
    assert 'undefined symbol' not in libraries_text and 'not found' not in libraries_text
    libraries = {str(Path(v).resolve(strict=True)) for v in re.findall(r'(/[^\s()]+)', libraries_text) if Path(v).is_file()}
    record['runtime_libraries'] = {p: pin(Path(p)) for p in sorted(libraries)}
    unchanged(prepared)
    record.update(passed=True, executable=identity(executable), link_proof=identity(W / 'link-proof.json'),
                  rename_proof=identity(W / 'rename-proof.json'), source_inputs_unchanged=True,
                  association_scope='Every retained function in renamed OLD Transport and exact new R Transport only. Other7 R native objects are pinned/linked, not falsely claimed covered by this two-object opcode proof.')
    write(result_path, record)
else:
    result_path = W / 'result.json'; assert not result_path.exists()
    build = json.loads((W / 'compile-result.json').read_text())
    assert build['passed'] and build['prepared'] == pin(W / 'prepared.json')
    for key in ('renamed_baseline_object', 'driver_object', 'executable', 'link_proof', 'rename_proof'):
        item = build[key]; assert pin(R / item['path']) == {k: item[k] for k in ('bytes', 'sha256')}
    for name, expected in build['runtime_libraries'].items():
        assert pin(Path(name)) == expected
    with (folder / 'comparison.json').open('x') as stream:
        code = subprocess.run([str(executable), str(folder)], cwd=R, stdout=stream).returncode
    try:
        comparison = json.loads((folder / 'comparison.json').read_text())
    except json.JSONDecodeError:
        comparison = {'passed': False, 'failure': 'missing_or_invalid_comparison'}
    hashes = {name: pin(folder / (name + '.bin')) for name in ('baseline', 'candidate') if (folder / (name + '.bin')).is_file()}
    before, after = folder / 'input-before.bin', folder / 'input-after.bin'
    input_equal = before.is_file() and after.is_file() and pin(before) == pin(after)
    first_input_mismatch = None
    if before.is_file() and after.is_file() and not input_equal:
        a, b = before.read_bytes(), after.read_bytes()
        offset = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
        first_input_mismatch = {'byte': offset, 'word': offset // 4,
                                'before_bytes': len(a), 'after_bytes': len(b)}
    cases = []
    if before.is_file() and (folder / 'cases.tsv').is_file():
        data = before.read_bytes()
        for line in (folder / 'cases.tsv').read_text().splitlines():
            index, family, name = line.split('\t'); index = int(index)
            row = data[(1 + index * 230) * 4:(1 + (index + 1) * 230) * 4]
            assert len(row) == 920
            cases.append({'index': index, 'family': family, 'name': name, 'bytes': len(row),
                          'sha256': hashlib.sha256(row).hexdigest()})
        assert len(data) == 100672 and struct.unpack_from('<I', data)[0] == 40
        program = data[(1 + 40 * 230) * 4:]
        plan = json.loads((R / 'out/native-cpu-accumulate-inline/coverage-plan.json').read_text())
        assert hashlib.sha256(program).hexdigest() == plan['inputs']['program']['sha256']
    counts = (comparison.get('case_count') == 38 and len(cases) == 38 and comparison.get('active_rows') == 38 and
              comparison.get('capacity_rows') == 40 and comparison.get('packets_per_implementation') == 5 and
              comparison.get('raw_calls') == 10 and comparison.get('compared_output_words') == 187960)
    complete_hashes = set(hashes) == {'baseline', 'candidate'} and all(v['bytes'] == 751840 for v in hashes.values())
    passed = (code == 0 and comparison.get('passed') and counts and input_equal and complete_hashes and
              hashes['baseline'] == hashes['candidate'] and comparison.get('input_unchanged') and
              comparison.get('inactive_sentinels_preserved') and comparison.get('full_word_match'))
    unchanged(prepared)
    record = {'passed': bool(passed), 'runtime_exit': code, 'comparison': comparison, 'cases': cases,
              'candidate_transport_object': prepared['candidate_transport_object'], 'all8_native_objects': prepared['objects'],
              'original_baseline_object': prepared['original_baseline_object'],
              'renamed_baseline_object': build['renamed_baseline_object'], 'executable': identity(executable),
              'link_proof': build['link_proof'], 'compile_result': identity(W / 'compile-result.json'),
              'output_hashes': hashes, 'input_before': pin(before) if before.is_file() else None,
              'input_after': pin(after) if after.is_file() else None, 'all_prepared_inputs_unchanged': True,
              'first_input_mismatch': first_input_mismatch,
              'scope': 'One finite compiler-change compatibility comparison of exact OLD/new raw Transport code, all38 frozen original inputs/all187960 words. No independent math, all-domain invalid-input, timing/benefit/frame/scientific/release/adoption claim. Original C8/ABBA22+Decimal/all8/frame remain mandatory.'}
    write(result_path, record)
    if not passed:
        raise SystemExit(1)
print(json.dumps({'mode': mode, 'passed': True}))
