from pathlib import Path
import hashlib, json, re, shutil, struct, subprocess

R = Path.cwd()
W = R / 'out/vulkan-descriptor-reuse/benefit'
T = R / 'bin/linux-gcc/vulkan-descriptor-reuse'
assert not W.exists()
W.mkdir()
P = R / 'attestations/software-vulkan/6014f91/portable-indexed-limb-accumulation/files/out/portable-indexed-limb-accumulation/benefit'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
write = lambda p, x: p.write_text(json.dumps(x, indent=2) + '\n')
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert subprocess.check_output(['git', 'status', '--porcelain'], text=True) == ''
baseline_revision = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
for name in ['software_contracts.py', 'compare_software.py', 'run_sequence.py']:
    shutil.copyfile(P / name, W / name)
for name in ['run_owned.py', 'owned_guard.py']:
    shutil.copyfile(W.parent / name, W / name)
wrapper = (P / 'run_limb_benefit_slot.py').read_text()
old = "[root/('bin/linux-gcc/portable-indexed-limb-accumulation/original_'+mode+'_controls')\n                      for mode in ('baseline','candidate')]"
new = "[root/'bin/linux-gcc/tests/backend/sirius_backend_tests',\n                      root/'bin/linux-gcc/vulkan-descriptor-reuse/original_candidate_controls']"
assert old in wrapper
(W / 'run_limb_benefit_slot.py').write_text(wrapper.replace(old, new))
arrays = {}
header = R / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
pattern = r'inline constexpr std::array<std::uint32_t,\s*(\d+)>\s+(\w+)\{\{(.*?)\}\};'
for size, name, body in re.findall(pattern, header.read_text(), re.S):
    words = [int(word) for word in re.findall(r'\b(\d+)u\b', body)]
    assert len(words) == int(size)
    data = struct.pack('<' + 'I' * len(words), *words)
    arrays[name] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
old_arrays = json.loads((P / 'baseline.json').read_text())['arrays']
assert len(arrays) == 40
assert all({k: old_arrays[name][k] for k in ('bytes', 'sha256')} == values for name, values in arrays.items())
for name in ['baseline.json', 'build-and-arrays.json']:
    write(W / name, {'arrays': arrays, 'compiled_header': {'path': str(header.relative_to(R)),
        'bytes': header.stat().st_size, 'sha256': sha(header)},
        'scope': 'All 40 original stock arrays extracted byte-exact from the actual unchanged primary GCC header. Both private producers use these exact modules; no shader overlay.'})
executables = {'baseline': R / 'bin/linux-gcc/tests/backend/sirius_backend_tests',
               'candidate': T / 'original_candidate_controls'}
before = json.loads((W.parent / 'preparation.json').read_text())['baseline_executable']
assert before['sha256'] == sha(executables['baseline']) and before['bytes'] == executables['baseline'].stat().st_size
old_actions = json.loads((P / 'actions.json').read_text())
base = ['/usr/lib/x86_64-linux-gnu/libvulkan.so', '/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so',
        '/usr/share/vulkan/icd.d/lvp_icd.json']
actions = {}
for label in ['a1', 'b1', 'b2', 'a2']:
    producer = 'baseline' if label.startswith('a') else 'candidate'
    resource_root = executables[producer].parent / 'resources'
    frozen = [str(p) for p in sorted(resource_root.rglob('*')) if p.is_file()]
    frozen += base + [str(executables[producer]), str(W / 'run_limb_benefit_slot.py'),
        str(W / 'software_contracts.py'), str(W / 'compare_software.py'), str(W / 'run_sequence.py'),
        str(W / 'baseline.json'), str(W / 'build-and-arrays.json'), str(W / 'retention-gate-frozen.json')]
    action = old_actions[label]
    actions[label] = {'revision': revision, 'kind': 'test', 'expected_exit': 0,
        'argv': ['python3', '-B', str(W / 'run_limb_benefit_slot.py'),
                 str(executables[producer].relative_to(R)), action['argv'][-1]],
        'limits': action['limits'], 'software_provider': True, 'frozen_inputs': frozen,
        'scope': 'Original unchanged 20+2 ABBA retained comparison, baseline exact accepted 4bcd987 ELF versus matching private descriptor candidate. No frame, full science or release claim.'}
write(W / 'actions.json', actions)
freeze = json.loads((P / 'retention-gate-frozen.json').read_text())
freeze.update(candidate_revision=revision, baseline_revision=baseline_revision,
    executables={k: {'path': str(p.relative_to(R)), 'source_revision': baseline_revision if k == 'baseline' else revision,
                     'bytes': p.stat().st_size, 'sha256': sha(p)} for k, p in executables.items()},
    actions_sha256=sha(W / 'actions.json'), scope='Same immutable reader, comparator, sequence and 20+2 criteria. Only executable path allowlist and actual bindings adapted; all numerical witnesses, repetitions and thresholds unchanged.')
assert freeze['reader_sha256'] == sha(W / 'software_contracts.py')
assert freeze['comparator_sha256'] == sha(W / 'compare_software.py')
assert freeze['launcher_sha256'] == sha(W / 'run_sequence.py')
write(W / 'retention-gate-frozen.json', freeze)
write(W / 'wrapper-adaptation.json', {'original_wrapper': {'path': str((P / 'run_limb_benefit_slot.py').relative_to(R)), 'sha256': sha(P / 'run_limb_benefit_slot.py')},
    'candidate_wrapper': {'path': str((W / 'run_limb_benefit_slot.py').relative_to(R)), 'sha256': sha(W / 'run_limb_benefit_slot.py')},
    'old_allowlist': old, 'new_allowlist': new, 'otherwise_byte_identical': True})
print(json.dumps({'sequence': freeze['sequence'], 'arrays': len(arrays), 'criteria_unchanged': True, 'executables': freeze['executables']}))
