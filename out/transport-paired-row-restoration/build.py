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
controller_signals=[]
def interrupted(signum,frame):controller_signals.append(signum)
signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
results = []
for label, args in commands:
    assert not (work / (label + '.stdout')).exists()
    started = time.monotonic()
    with (work / (label + '.stdout')).open('wb') as out, (work / (label + '.stderr')).open('wb') as err:
        assert not controller_signals
        child = subprocess.Popen(args, stdout=out, stderr=err, start_new_session=True)
        try:
            birth = Path('/proc', str(child.pid), 'stat').read_text().rsplit(')',1)[1].split()[19]
            while child.poll() is None:
                assert not controller_signals, controller_signals
                remaining=600-(time.monotonic()-started)
                if remaining<=0:raise subprocess.TimeoutExpired(args,600)
                try:child.wait(timeout=min(0.25,remaining))
                except subprocess.TimeoutExpired:pass
            code = child.returncode
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
    assert code == 0 and not controller_signals, label
header = root / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
arrays = {}
for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', header.read_text(), re.S):
    words = [int(q[:-1], 0) for q in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u', body)]
    assert len(words) == int(count)
    payload = struct.pack('<'+'I'*len(words), *words)
    arrays[name] = dict(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
assert len(arrays) == 40
changed = [name for name in arrays if arrays[name] != baseline['arrays'][name]]
assert changed == [], changed
assert arrays == baseline['arrays']
header_sha = hashlib.sha256(header.read_bytes()).hexdigest()
assert header_sha == baseline['whole_header_sha256']
protected_header = root / 'attestations/native-vulkan/6845c81/paired-transport-source-numerical/source-inputs/bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
assert hashlib.sha256(protected_header.read_bytes()).hexdigest() == baseline['whole_header_sha256']
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == head
assert not subprocess.check_output(['git', 'status', '--porcelain'])
record = dict(kind='Restored accepted Transport after rejected paired-row gate', source_revision=head,
              source_tree=subprocess.check_output(['git', 'rev-parse', 'HEAD^{tree}'], text=True).strip(),
              commands=results, returncode=0, restored_all_40_arrays_exact=True,
              whole_header_sha256=header_sha, arrays=arrays,
              scope='Strict existing preset source/alignment/app/render/backend restoration build; whole generated header and all40 complete arrays exact accepted baseline. No numerical execution/frame/fullqualification claim.')
(work / 'build-and-arrays.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(dict(source_revision=head, restored_all_40_arrays_exact=True, header_sha256=header_sha,
                      elapsed_seconds=sum(q['elapsed_seconds'] for q in results))))
