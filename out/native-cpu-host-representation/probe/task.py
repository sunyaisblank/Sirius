"""Exactly one compiler/link invocation or one generated comparison executable."""
from pathlib import Path
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2); stream.write('\n')


mode, output = sys.argv[1:]
assert not any(os.environ.get(key) for key in
               ('CPATH', 'CPLUS_INCLUDE_PATH', 'C_INCLUDE_PATH', 'LIBRARY_PATH', 'GCC_EXEC_PREFIX',
                'COMPILER_PATH', 'GCC_SPECS', 'LD_PRELOAD', 'LD_AUDIT', 'LD_LIBRARY_PATH'))
folder = Path(output); assert folder == W / mode and mode in ('compile', 'run')
prepared = json.loads((W / 'preparation.json').read_text())
for name, value in prepared['tool_dependencies'].items():
    assert pin(Path(name)) == value
for name, value in prepared['direct_headers'].items():
    assert pin(Path(name)) == value
generation = json.loads((W / 'generation.json').read_text())
assert pin(W / 'probe.cpp') == generation['generated_source']
for key in ('source', 'patch'):
    assert pin(R / generation[key]) == generation[key + '_identity']
executable = R / prepared['executable']
if mode == 'compile':
    assert not executable.exists() and not (W / 'compile-result.json').exists()
    with (folder / 'compiler.log').open('x') as log:
        code = subprocess.run(prepared['argv'], cwd=R, stdout=log, stderr=subprocess.STDOUT).returncode
    receipt = {'passed': False, 'exit': code, 'argv': prepared['argv'], 'compiler_log': pin(folder / 'compiler.log'),
               'generation': pin(W / 'generation.json'), 'preparation': pin(W / 'preparation.json')}
    if code == 0:
        dependency_text = (folder / 'dependencies.d').read_text().replace('\\\n', ' ')
        paths = shlex.split(dependency_text.split(':', 1)[1])
        headers = {str((Path(path) if Path(path).is_absolute() else R / path).resolve(strict=True)) for path in paths}
        receipt['headers'] = {path: pin(Path(path)) for path in sorted(headers)}
        text = subprocess.check_output([prepared['ldd'], str(executable)], text=True, timeout=10)
        (folder / 'linked-libraries.txt').write_text(text)
        libraries = {str(Path(match).resolve(strict=True)) for match in re.findall(r'(/[^\s]+)\s+\(', text)}
        receipt['runtime_libraries'] = {path: pin(Path(path)) for path in sorted(libraries)}
        for path, value in receipt['runtime_libraries'].items():
            assert prepared['tool_dependencies'].get(path) == value, ('unexpected runtime dependency', path)
        linked = set()
        for line in (folder / 'compiler.log').read_text().splitlines():
            path = Path(line.strip())
            if path.is_absolute() and path.is_file():
                linked.add(str(path.resolve(strict=True)))
        receipt.update(executable=pin(executable), link_inputs={p: pin(Path(p)) for p in sorted(linked)}, passed=True)
    write(W / 'compile-result.json', receipt)
    raise SystemExit(0 if receipt['passed'] else 1)

compiled = json.loads((W / 'compile-result.json').read_text())
assert compiled['passed'] and pin(executable) == compiled['executable']
for category in ('headers', 'runtime_libraries', 'link_inputs'):
    for name, value in compiled[category].items():
        assert pin(Path(name)) == value
assert not (W / 'result.json').exists()
with (folder / 'comparison.json').open('x') as log, (folder / 'probe-stderr.log').open('x') as errors:
    code = subprocess.run([str(executable), str(folder / 'actual-maps.txt')], cwd=R, stdout=log, stderr=errors).returncode
raw = (folder / 'comparison.json').read_text()
try:
    comparison = json.loads(raw)
except json.JSONDecodeError:
    comparison = {'passed': False, 'invalid_or_missing_comparison': True}
mapped = set()
if (folder / 'actual-maps.txt').exists():
    for line in (folder / 'actual-maps.txt').read_text().splitlines():
        fields = line.split(None, 5)
        if len(fields) == 6 and fields[5].startswith('/'):
            mapped.add(str(Path(fields[5]).resolve(strict=True)))
identities = {p: pin(Path(p)) for p in sorted(mapped)}
expected = {**compiled['runtime_libraries'], str(executable): compiled['executable']}
counts = (comparison.get('unique_inputs') == prepared['expected_unique_inputs'] and
          comparison.get('environments') == prepared['expected_environments'] and
          comparison.get('comparisons') == prepared['expected_comparisons'] and
          comparison.get('function_calls') == prepared['expected_comparisons'] * 2)
passed = code == 0 and comparison.get('passed') and comparison.get('ambient_restored') and counts and identities == expected
receipt = {'passed': bool(passed), 'exit': code, 'comparison': comparison, 'exact_expected_counts': counts,
           'actual_mapped_files': identities, 'all_actual_runtime_bytes_match_compile_pins': identities == expected,
           'executable': pin(executable), 'generated_source': pin(W / 'probe.cpp'),
           'generation': pin(W / 'generation.json'), 'preparation': pin(W / 'preparation.json'),
           'raw_comparison': pin(folder / 'comparison.json'), 'stderr': pin(folder / 'probe-stderr.log'),
           'scope': 'This exact binary/compiler/libm, masked-trap32-environment finite boundary matrix only. Not exhaustive binary32, all-platform fenv/errno or performance/scientific qualification.'}
write(W / 'result.json', receipt)
print(json.dumps({'passed': receipt['passed'], 'comparisons': comparison.get('comparisons'), 'exit': code}))
raise SystemExit(0 if receipt['passed'] else 1)
