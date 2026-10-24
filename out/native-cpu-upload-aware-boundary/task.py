"""Execute exactly one prepared fixture build or one fake-device control."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


def unchanged(prepared):
    for name, expected in prepared['source_inputs'].items():
        assert pin(R / name) == expected, name


mode, output = sys.argv[1:]
assert mode in ('compile', 'run') and Path(output) == W / mode
assert not any(os.environ.get(k) for k in ('CPATH', 'CPLUS_INCLUDE_PATH', 'LIBRARY_PATH', 'GCC_EXEC_PREFIX', 'COMPILER_PATH'))
prepared = json.loads((W / 'prepared.json').read_text())
unchanged(prepared)
if mode == 'compile':
    destination = W / 'compile-result.json'
    assert not destination.exists()
    result = {'passed': False, 'commands': [], 'extra_fixture_CPP_compilations': 0,
              'extra_fixture_links': 0, 'observer_compilations': 0, 'P_host_rebuilds': 0,
              'shader_provider_compilations': 0, 'prepared': pin(W / 'prepared.json')}
    for index, command in enumerate(prepared['commands']):
        code = subprocess.run(command['argv'], cwd=command['cwd']).returncode
        result['commands'].append({**command, 'exit': code})
        result['extra_fixture_CPP_compilations' if index == 0 else 'extra_fixture_links'] += 1
        destination.write_text(json.dumps(result, indent=2) + '\n')
        if code:
            raise SystemExit(code)
    unchanged(prepared)
    for key in ('object', 'executable'):
        result[key] = {'path': prepared[key], **pin(R / prepared[key])}
    result.update(passed=True, backend=prepared['backend'], backend_identity=prepared['backend_identity'],
                  all_prepared_source_inputs_unchanged=True)
    destination.write_text(json.dumps(result, indent=2) + '\n')
else:
    destination = W / 'result.json'
    assert not destination.exists()
    built = json.loads((W / 'compile-result.json').read_text())
    assert built['passed'] and built['prepared'] == pin(W / 'prepared.json')
    assert pin(R / prepared['executable']) == {k: built['executable'][k] for k in ('bytes', 'sha256')}
    raw = Path(output) / 'wire-result.json'
    with raw.open('x') as stream:
        code = subprocess.run([str(R / prepared['executable'])], cwd=R, stdout=stream).returncode
    try:
        comparison = json.loads(raw.read_text())
    except json.JSONDecodeError:
        comparison = {'passed': False, 'invalid_or_missing_result': True}
    unchanged(prepared)
    passed = (code == 0 and comparison.get('passed') and comparison.get('write_attempts') == 4 and
              comparison.get('successful_dispatches') == comparison.get('successful_reads') == 3 and
              comparison.get('first_success_zero_inactive') and comparison.get('later_prefix_preserved_inactive'))
    result = {'passed': bool(passed), 'exit': code, 'wire_result': comparison, 'raw_result': pin(raw),
              'fixture_executions': 1, 'compile_result': pin(W / 'compile-result.json'),
              'prepared': pin(W / 'prepared.json'), 'executable': pin(R / prepared['executable']),
              'backend': prepared['backend'], 'backend_identity': prepared['backend_identity'],
              'all_prepared_source_inputs_unchanged': True,
              'scope': 'Only the actual P backend failure/publication/prefix-upload contract against one fake device. No physical correctness, frame/performance, all8/all22, scientific/release or adoption claim.'}
    with destination.open('x') as stream:
        json.dump(result, stream, indent=2); stream.write('\n')
    if not passed:
        raise SystemExit(1)
print(json.dumps({'mode': mode, 'passed': True}))
