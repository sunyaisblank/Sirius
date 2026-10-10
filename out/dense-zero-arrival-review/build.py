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
source = json.loads((work/'source.json').read_text())
assert head == source['source_revision']
restored = source['trial_disposition'] == 'rejected_restored'
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
changed = {name for name in arrays if arrays[name] != baseline['arrays'][name]}
assert changed == (set() if restored else {'kDenseShader', 'kDenseFp64Shader', 'kDensePortableShader', 'kDensePortableFp64Shader'}), changed
header_sha = hashlib.sha256(header.read_bytes()).hexdigest()
assert (header_sha == baseline['whole_header_sha256']) == restored
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == head
assert not subprocess.check_output(['git', 'status', '--porcelain'])
record = dict(kind='Restored Dense source with retained zero-arrival regression' if restored else 'Canonical exact-zero moving Dense arrival trial', source_revision=head,
              source_tree=subprocess.check_output(['git', 'rev-parse', 'HEAD^{tree}'], text=True).strip(),
              commands=results, returncode=0, all_36_unaffected_arrays_exact=True, all_40_arrays_exact=restored, changed_array_names=sorted(changed),
              whole_header_sha256=header_sha, arrays=arrays,
              scope=('Clean exact-source strict source/alignment/app/render/backend build. All40 arrays and whole header are exactly restored to the accepted baseline; useful independent regression and Mandatory governance retained. No device/fullscience/release claim.' if restored else 'Clean exact-source strict source/alignment/app/render/backend build of the Dense arrival trial. All36 unaffected arrays and all original program tables remain exact; exactly four Dense shader modules change. No device, benefit, fullscience or release claim.'))
(work / 'build-and-arrays.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(dict(source_revision=head, all_36_unaffected_arrays_exact=True, all_40_arrays_exact=restored, changed_array_names=sorted(changed), header_sha256=header_sha,
                      elapsed_seconds=sum(q['elapsed_seconds'] for q in results))))
