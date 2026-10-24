"""Seal existing compiled products and the sole C8; never execute it here."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
native = root / 'out/native-cpu-slang-bridge'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


assert not (work / 'actions.json').exists()
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
reuse = doc(native / 'reused-row-products.json')
build = doc(native / 'provider/provider-build-result.json')
compiled = doc(native / 'batch-compile-result.json')
assert build['passed'] and doc(native / 'provider/linked-native-proof.json')['passed']
assert doc(native / 'ordinary/guard-acceptance.json')['passed']
executable = root / 'bin/linux-gcc/native-cpu-row-ownership/matched/audited_backend_tests'
expected = {name: {k: value[k] for k in ('bytes', 'sha256')}
            for name, value in reuse['restored_actual_products'].items()}
expected.update(compiled['source_inputs'])
expected.update({item['path']: {k: item[k] for k in ('bytes', 'sha256')}
                 for item in build['native_objects']})
library, icd = root / build['library'], root / build['ICD']
expected[str(library)] = build['library_identity']
expected[str(icd)] = build['ICD_identity']
for name, value in expected.items():
    assert pin(root / name) == value, name
assert pin(executable)['sha256'] == '7f2d1fa88953b954ab116233458a4cb962470f863fda9bbb3a1cad7d4a7d6a26'
files = set(expected)
files |= {str(p.relative_to(root)) for p in work.iterdir() if p.is_file()}
files |= {str(p.relative_to(root)) for p in (root / 'bin/linux-gcc/src/sirius/backend/retained').iterdir() if p.is_file()}
files |= {str(p.relative_to(root)) for p in (native / 'provider').iterdir() if p.is_file()}
files |= {'out/native-cpu-row-ownership/compiler-result.json',
          'out/native-cpu-row-ownership/cold-frame/run_owned.py',
          'out/native-cpu-row-ownership/cold-frame/owned_guard.py',
          'out/native-cpu-arithmetic-sharing/readers/native_evidence.py',
          'out/native-cpu-arithmetic-sharing/readers/expected-readback-layout.json',
          'out/native-cpu-row-ownership/previous-control-8-acceptance.json',
          'out/native-original-input-ownership/provider-history/candidate-default/stdout.log',
          'out/native-cpu-slang-bridge/ordinary/transport/device-features.json',
          'out/native-cpu-slang-bridge/ordinary/transport/float-controls-properties.json',
          '/etc/ld.so.cache'}
for path in [executable, library, Path('/usr/bin/gdb'), Path('/usr/bin/python3'), Path('/usr/bin/nm')]:
    files |= {str(path), str(path.resolve(strict=True))}
    output = subprocess.check_output(['/usr/bin/ldd', str(path)], text=True)
    files |= {str(Path(name).resolve(strict=True)) for name in re.findall(r'(/[^\s()]+)', output)
              if Path(name).is_file()}
assert all((root / name).is_file() for name in files)
action = {'arm': 'candidate', 'case': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults',
    'prerequisites': [], 'kind': 'test', 'expected_exit': 0,
    'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'argv': ['/usr/bin/python3', '-B', str(work / 'run_case.py'), 'control-8', '{folder}'],
    'frozen_inputs': sorted(files), 'expected_identities': expected,
    'executable': str(executable), 'software_provider': True,
    'runtime_icd': str(icd), 'runtime_provider': str(library),
    'scope': 'One original C8 only, unchanged complete-row matched ELF7f2d1fa8 and nativeCPU providerccb9115a. First actual native packet/return perstage observed passively; all original counts/assertions/cache/governors/300s4GiB2GiB unchanged. C8 timing under debugger is not benefit evidence. All32complete digests/51markers/54modulejoins must pass original reader. No ABBA/frame/qualification/adoption follows unless separately admitted.'}
(work / 'actions.json').write_text(json.dumps({'control-8': action}, indent=2) + '\n')
print(json.dumps({'source': action['revision'], 'frozen_inputs': len(files), 'executed': False}))
