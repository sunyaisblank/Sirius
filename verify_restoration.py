"""Verify the terminal restored build against preserved complete original payloads."""
from pathlib import Path
import hashlib, json, re, struct, subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
checks = work / 'restoration'
read = lambda p: json.loads(p.read_text())
head = read(checks / 'source.json')['source_revision']
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == head
assert not subprocess.check_output(['git', 'status', '--porcelain'])
assert subprocess.check_output(['git', 'rev-parse', 'HEAD^{tree}']) == subprocess.check_output(['git', 'rev-parse', '7f028007a4935ffdefe5ee9c886b138299f50dd8^{tree}'])
for case in ('configure', 'build'):
    owner = read(checks / case / 'owner.json')
    assert owner['passed'] and owner['child_exit'] == 0
    assert owner['source'] == {'revision': head, 'status': ''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and not owner['cleanup_errors']
    assert not owner['remaining_owned_processes'] and all(v['absent'] for v in owner['observed_births'])
    assert read(checks / case / 'inputs-before.json') == read(checks / case / 'inputs-after.json')
header = root / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
data = header.read_bytes()
arrays = {}
for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', data.decode(), re.S):
    words = [int(q[:-1], 0) for q in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u', body)]
    assert name not in arrays and len(words) == int(count)
    payload = struct.pack('<' + 'I' * len(words), *words)
    arrays[name] = {'bytes': len(payload), 'sha256': hashlib.sha256(payload).hexdigest()}
baseline = read(work / 'baseline.json')
old_restoration = read(root / 'attestations/review/7f02800/rejected-fma-restoration-source-checks/diagnostic/restoration/build-and-arrays.json')
sha = hashlib.sha256(data).hexdigest()
assert len(arrays) == 40 and arrays == baseline['arrays'] == old_restoration['arrays']
assert sha == baseline['whole_header_sha256'] == old_restoration['whole_header_sha256']
assert old_restoration['original40_complete_arrays_and_whole_header_byte_exact95']
expected_spvs = {root / p for p in read(checks / 'actions.json')['admission']['frozen_inputs'] if p.endswith('.spv')}
assert len(expected_spvs) == 33 and set(header.parent.glob('*.spv')) == expected_spvs
record = {'source_revision': head, 'entire_tree_byte_exact_7f': True,
          'whole_header_sha256': sha, 'arrays': arrays,
          'original40_complete_arrays_and_whole_header_byte_exact95': True,
          'declared_spv_count': 33,
          'scope': 'Terminal strict source/build checks and complete original header/arrays only. Actual consumer joins and two host admission cases follow separately; no current Vulkan, frame, science or full qualification execution.'}
(checks / 'build-and-arrays.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'source_revision': head, 'complete_arrays': 40, 'whole_header_sha256': sha, 'spvs': 33}))
