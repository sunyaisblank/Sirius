from pathlib import Path
import hashlib
import json
import re
import struct
import subprocess
import time

root = Path.cwd()
work = Path(__file__).resolve().parent
head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert head == json.loads((work/'candidate-source.json').read_text())['source_revision']
assert not subprocess.check_output(['git', 'status', '--porcelain'])
baseline = json.loads((work/'baseline-build-and-arrays.json').read_text())
plan = json.loads((work/'candidate-plan-binding.json').read_text())
assert hashlib.sha256((root/'src/sirius/kernels/retained_transport.slang').read_bytes()).hexdigest() == plan['candidate_source_sha256']
assert hashlib.sha256((root/'scripts/build-retained-kernels.py').read_bytes()).hexdigest() == plan['candidate_builder_sha256']
commands = [
    ('candidate-configure', ['cmake', '--preset', 'linux-gcc', '-DSIRIUS_ALIGNMENT_MODE=qualification']),
    ('candidate-build', ['cmake', '--build', '--preset', 'linux-gcc', '--target', 'SiriusAlignmentGate', 'SiriusSourceGovernance', 'sirius', 'sirius_app_tests', 'sirius_render_tests', 'sirius_backend_tests', '-j2']),
]
results = []
for label, args in commands:
    assert not (work / (label + '.stdout')).exists()
    started = time.monotonic()
    with (work / (label + '.stdout')).open('wb') as out, (work / (label + '.stderr')).open('wb') as err:
        child = subprocess.Popen(args, stdout=out, stderr=err)
        birth = Path('/proc', str(child.pid), 'stat').read_text().split(') ', 1)[1].split()[19]
        code = child.wait()
    results.append(dict(label=label, arguments=args, pid=child.pid, start_ticks=birth,
                        returncode=code, elapsed_seconds=time.monotonic()-started))
    (work / 'candidate-build-commands.json').write_text(json.dumps(results, indent=2) + '\n')
    assert code == 0, label
header = root / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
arrays = {}
for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', header.read_text(), re.S):
    words = [int(q[:-1], 0) for q in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u', body)]
    assert len(words) == int(count)
    payload = struct.pack('<'+'I'*len(words), *words)
    arrays[name] = dict(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
assert len(arrays) == 40 and set(arrays)==set(baseline['arrays'])
changed = sorted(name for name in arrays if arrays[name] != baseline['arrays'][name])
assert changed == sorted(plan['changed_expected_arrays']), changed
assert len(arrays)-len(changed)==34
header_sha = hashlib.sha256(header.read_bytes()).hexdigest()
assert header_sha != baseline['whole_header_sha256']
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == head
assert not subprocess.check_output(['git', 'status', '--porcelain'])
record = dict(kind='Homogeneous Transport opcode-dispatch bounded candidate build', source_revision=head,
              source_tree=subprocess.check_output(['git', 'rev-parse', 'HEAD^{tree}'], text=True).strip(),
              commands=results, returncode=0, all_40_arrays_bound=True, excluded_34_arrays_exact=True, changed_arrays=changed,
              whole_header_sha256=header_sha, arrays=arrays,
              scope='Existing source/alignment/app/render/backend candidate build; exactly six Transport shader arrays changed, all other34/table bytes exact. No numerical execution, performance or full qualification claim.')
(work / 'candidate-build-and-arrays.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(dict(source_revision=head, all_40_arrays_bound=True, excluded_34_arrays_exact=True, changed_arrays=changed, header_sha256=header_sha,
                      elapsed_seconds=sum(q['elapsed_seconds'] for q in results))))
