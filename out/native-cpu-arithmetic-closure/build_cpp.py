"""Build the sealed C++ helper/reference fixture once, preserving commands."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
prepared = json.loads((W / 'prepared-cpp.json').read_text())


def pin(path):
    return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}


def verify():
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == prepared['revision']
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    for name, expected in prepared['inputs'].items():
        assert pin(R / name) == expected, name


verify()
destination = (R / prepared['executable']).parent
assert not destination.exists() and not (W / 'cpp-build-result.json').exists()
destination.mkdir(parents=True)
result = {'passed': False, 'compilations': [], 'new_CPP_compilations': 0, 'new_links': 0}


def save():
    (W / 'cpp-build-result.json').write_text(json.dumps(result, indent=2) + '\n')


save()
for item in prepared['compilations']:
    row = {'name': item['name'], 'argv': item['argv'], 'started': True}
    result['compilations'].append(row)
    save()
    code = subprocess.run(item['argv'], cwd=R).returncode
    row['exit'] = code
    result['new_CPP_compilations'] += 1
    save()
    if code:
        raise SystemExit(code)
    row['object'] = {'path': item['object'], **pin(R / item['object'])}
    save()
code = subprocess.run(prepared['link'], cwd=R).returncode
result['new_links'] = 1
result['link_argv'], result['link_exit'] = prepared['link'], code
save()
if code:
    raise SystemExit(code)
executable = R / prepared['executable']
result['executable'] = {'path': prepared['executable'], **pin(executable)}
runtime = subprocess.check_output(['/usr/bin/ldd', '-r', str(executable)], text=True, stderr=subprocess.STDOUT)
assert 'undefined symbol' not in runtime and 'not found' not in runtime
result['runtime_libraries'] = {str(Path(name).resolve()): pin(Path(name).resolve())
    for name in re.findall(r'(/[^\s()]+)', runtime) if Path(name).is_file()}
verify()
result['passed'] = True
save()
print(json.dumps({'passed': True, 'new_CPP_compilations': 2, 'new_links': 1}))
