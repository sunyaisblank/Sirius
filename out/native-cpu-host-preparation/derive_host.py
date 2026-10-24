"""Derive one private host copy with exact inverse to the matched7f2 source."""
from pathlib import Path
import hashlib
import json
import re
import sys

root = Path.cwd()
work = Path(__file__).resolve().parent
source_path = root / 'out/native-cpu-row-ownership/retained_compute.cpp'
destination = root / 'bin/linux-gcc/native-cpu-host-preparation/source/retained_compute.cpp'
assert sys.argv[1:] in ([], ['--write'])


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        assert sys.argv[1:] == ['--write'], path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)


copy_proof_path = root / 'out/native-cpu-host-copy-prepare/proof.json'
copy_proof = json.loads(copy_proof_path.read_text())
source = source_path.read_text()
assert pin(source.encode()) == copy_proof['sources'][str(source_path.relative_to(root))]['basis']
spacing_patch = root / 'out/native-cpu-host-representation/host.patch'
lines = spacing_patch.read_text().splitlines(True)
assert sum(line.startswith('@@') for line in lines) == 1
hunk = lines[next(i for i, line in enumerate(lines) if line.startswith('@@')) + 1:]
assert all(line[:1] in (' ', '+', '-') for line in hunk)
old_spacing = ''.join(line[1:] for line in hunk if line[:1] in (' ', '-'))
new_spacing = ''.join(line[1:] for line in hunk if line[:1] in (' ', '+'))
assert source.count(old_spacing) == 1
text = source.replace(old_spacing, new_spacing, 1)
fallback = ('return double(std::nextafter(magnitude, std::numeric_limits<float>::infinity())) -\n'
            '               double(magnitude);')
assert old_spacing.count(fallback) == new_spacing.count(fallback) == 1
edits, sites = [], []
pattern = re.compile(r'    std::fill_n\((\w+_)\.input.begin\(\) \+ 1, capacity_ \* (\d+), 0U\);\n'
                     r'    std::memcpy\(\1\.input.data\(\) \+ 1, (\w+)\.data\(\), \3\.size_bytes\(\)\);')
for match in pattern.finditer(text):
    stage, stride, span = match.groups()
    replacement = (f'    std::fill_n({stage}.input.begin() + 1 + {span}.size() * {stride},\n'
                   f'                (capacity_ - {span}.size()) * {stride}, 0U);\n'
                   f'    std::memcpy({stage}.input.data() + 1, {span}.data(), {span}.size_bytes());')
    edits.append((match[0], replacement))
    sites.append({'stage': stage, 'words_per_row': int(stride), 'input': span})
assert len(edits) == 9
assert sites == [{key: value for key, value in item.items() if key != 'line'}
                 for item in copy_proof['sources'][str(source_path.relative_to(root))]['sites']]
for old, new in edits:
    assert text.count(old) == 1
    text = text.replace(old, new, 1)
restored = text
for old, new in reversed(edits):
    assert restored.count(new) == 1
    restored = restored.replace(new, old, 1)
assert restored.count(new_spacing) == 1
restored = restored.replace(new_spacing, old_spacing, 1)
assert restored == source
inputs = {str(path.relative_to(root)): pin(path.read_bytes())
          for path in (source_path, spacing_patch, copy_proof_path,
                       root / 'out/native-cpu-host-copy-prepare/matched-host.patch')}
libm = Path('/usr/lib/x86_64-linux-gnu/libm.so.6')
inputs[str(libm)] = pin(libm.read_bytes())
proof = {'basis_path': str(source_path.relative_to(root)), 'basis': pin(source.encode()),
    'derived_path': str(destination.relative_to(root)), 'derived': pin(text.encode()),
    'source_inputs': inputs, 'full_inverse_restores_matched_source': True,
    'normal_spacing_patch_applied_once_in_private_copy': True,
    'original_spacing_fallback_verbatim': True, 'inactive_only_input_zero_ranges': sites,
    'all_other_source_bytes_unchanged': True, 'source_selectors_and_all_shader_arrays_unchanged': True,
    'normal_spacing_scope': 'Original abs remains first. Only finite positive normal magnitudes below maxfinite use the exact exponent-derived normal binary64 ULP. All zero/subnormal/maxfinite/nonfinite cases retain the original nextafter expression and its exception/errno behavior. Root independently reviews/tests linked-libm FEnv/errno equivalence before runtime admission.',
    'copy_scope': 'Complete active memcpy retained, inactive rows still zero, header/program suffix unchanged. Final host staging/device upload bytes identical for every validated size.',
    'limits': 'Derived source only; no compiler/build/runtime, bottleneck, speedup, scientific or release claim.'}
publish(destination, text.encode())
publish(work / 'host-source-proof.json', (json.dumps(proof, indent=2) + '\n').encode())
assert source_path.read_text() == source
print(json.dumps({'private_source_derived': True, 'zero_range_changes': 9,
                  'spacing_changes': 1, 'full_inverse_exact': True, 'compiled_or_executed': False}))
