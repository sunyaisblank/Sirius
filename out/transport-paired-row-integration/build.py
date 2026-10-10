from pathlib import Path
import hashlib
import json
import re
import struct
import subprocess
import time
import os
import signal

root = Path.cwd()
work = Path(__file__).resolve().parent
head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert head == json.loads((work/'source.json').read_text())['source_revision']
assert not subprocess.check_output(['git', 'status', '--porcelain'])
baseline = json.loads((root / 'attestations/native-vulkan/9d56406/endpoint-native-lanes32-rejected/baseline-arrays.json').read_text())
commands = [
    ('configure', ['cmake', '--preset', 'linux-gcc', '-DSIRIUS_ALIGNMENT_MODE=qualification']),
    ('build', ['cmake', '--build', '--preset', 'linux-gcc', '--target', 'SiriusAlignmentGate', 'SiriusSourceGovernance', 'sirius', 'sirius_app_tests', 'sirius_render_tests', 'sirius_backend_tests', '-j2']),
]
results = []
for label, args in commands:
    assert not (work / (label + '.stdout')).exists()
    started = time.monotonic()
    with (work / (label + '.stdout')).open('wb') as out, (work / (label + '.stderr')).open('wb') as err:
        child = subprocess.Popen(args, stdout=out, stderr=err, start_new_session=True)
        birth = Path('/proc', str(child.pid), 'stat').read_text().split(') ', 1)[1].split()[19]
        try:
            code = child.wait(timeout=600)
        finally:
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGTERM)
                try: child.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait(timeout=10)
    results.append(dict(label=label, arguments=args, pid=child.pid, start_ticks=birth,
                        returncode=code, elapsed_seconds=time.monotonic()-started))
    (work / 'build-commands.json').write_text(json.dumps(results, indent=2) + '\n')
    assert code == 0, label
header = root / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
arrays = {}
for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', header.read_text(), re.S):
    words = [int(q[:-1], 0) for q in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u', body)]
    assert len(words) == int(count)
    payload = struct.pack('<'+'I'*len(words), *words)
    arrays[name] = dict(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
assert len(arrays) == 40
changed = [name for name in arrays if arrays[name] != baseline['arrays'][name]]
assert changed == ['kTransportShader'], changed
assert arrays['kTransportShader'] == {'bytes': 92080, 'sha256': '721be7625f8eebbf722adf3e79a89050bffe9bf9b7260b8f5120b0a85ccd7f27'}
header_sha = hashlib.sha256(header.read_bytes()).hexdigest()
protected_header = root / 'attestations/native-vulkan/6845c81/paired-transport-source-numerical/source-inputs/bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
assert hashlib.sha256(protected_header.read_bytes()).hexdigest() == baseline['whole_header_sha256']
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == head
assert not subprocess.check_output(['git', 'status', '--porcelain'])
record = dict(kind='Unadopted native-default paired Transport integration', source_revision=head,
              source_tree=subprocess.check_output(['git', 'rev-parse', 'HEAD^{tree}'], text=True).strip(),
              commands=results, returncode=0, reviewed_native_transport_exact=True, other_39_arrays_exact=True,
              whole_header_sha256=header_sha, arrays=arrays,
              scope='Strict existing preset source/alignment/app/render/backend build at clean unadopted candidate; reviewed native Transport721 exact and other39 accepted arrays exact. No numerical execution, performance/adoption/full qualification claim.')
(work / 'build-and-arrays.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(dict(source_revision=head, reviewed_native_transport_exact=True, other_39_arrays_exact=True, header_sha256=header_sha,
                      elapsed_seconds=sum(q['elapsed_seconds'] for q in results))))
