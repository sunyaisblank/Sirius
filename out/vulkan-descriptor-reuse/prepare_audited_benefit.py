from pathlib import Path
import hashlib, json, shlex, shutil, subprocess

R = Path.cwd()
Q = Path(__file__).resolve().parent
W = Q / 'audited-benefit'
T = R / 'bin/linux-gcc/vulkan-descriptor-reuse'
P = R / 'attestations/software-vulkan/6014f91/portable-indexed-limb-accumulation/files/out/portable-indexed-limb-accumulation/benefit'
assert not W.exists()
W.mkdir()
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
write = lambda p, x: p.write_text(json.dumps(x, indent=2) + '\n')
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert revision == '8a27321c53332fb76273556987040e9ac6ab8fa2'
assert subprocess.check_output(['git', 'status', '--porcelain'], text=True) == ''
for name in ['software_contracts.py', 'compare_software.py', 'run_sequence.py', 'expected-readback-layout.json']:
    shutil.copyfile(P / name, W / name)
for name in ['run_owned.py', 'owned_guard.py']:
    shutil.copyfile(Q / name, W / name)
for name in ['baseline.json', 'build-and-arrays.json']:
    shutil.copyfile(Q / 'benefit' / name, W / name)
old = "[root/('bin/linux-gcc/portable-indexed-limb-accumulation/original_'+mode+'_controls')\n                      for mode in ('baseline','candidate')]"
new = "[root/('bin/linux-gcc/vulkan-descriptor-reuse/audited_'+mode+'_controls')\n                      for mode in ('baseline','candidate')]"
wrapper = (P / 'run_limb_benefit_slot.py').read_text()
assert old in wrapper
(W / 'run_limb_benefit_slot.py').write_text(wrapper.replace(old, new))
executables = {mode: T / ('audited_' + mode + '_controls') for mode in ['baseline', 'candidate']}
provider = ['/usr/lib/x86_64-linux-gnu/libvulkan.so', '/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so',
            '/usr/share/vulkan/icd.d/lvp_icd.json']
resources = [str(p) for p in sorted((T / 'resources').rglob('*')) if p.is_file()]
actions = {}
maps = {}
for mode in ['baseline', 'candidate']:
    link = json.loads((Q / ('link-audited-' + mode) / 'owner.json').read_text())
    assert link['passed'] and link['status'] == 'terminal'
    path = T / ('audited-' + mode + '.map')
    text = path.read_text()
    assert 'libsirius_backend.a(vulkan_device.cpp.o)' not in text
    assert 'vulkan_smoke_test.cpp.o' not in text
    assert str(T / ('audited-' + mode + '-observer.o')) in text
    assert 'audited-' + ('candidate' if mode == 'baseline' else 'baseline') + '-observer.o' not in text
    adapter = str(R / ('bin/linux-gcc/src/sirius/backend/CMakeFiles/sirius_backend.dir/vulkan/vulkan_device.cpp.o'
                       if mode == 'baseline' else 'bin/linux-gcc/vulkan-descriptor-reuse/vulkan_device.o'))
    assert adapter in text
    dep = (T / ('audited-' + mode + '-observer.d')).read_text()
    actual_headers = [str(Path(word).resolve()) for word in
                      shlex.split(dep.replace('\\\n', '').split(':', 1)[1])]
    header = str(T / 'baseline-header/sirius/backend/vulkan/vulkan_device.h') if mode == 'baseline' else str(R / 'src/sirius/backend/vulkan/vulkan_device.h')
    assert header in actual_headers
    maps[mode] = {'path': str(path.relative_to(R)), 'sha256': sha(path), 'adapter': adapter,
                  'concrete_header': header, 'single_explicit_observer_and_adapter': True}
for label in ['a1', 'b1', 'b2', 'a2']:
    mode = 'baseline' if label.startswith('a') else 'candidate'
    original = json.loads((P / 'actions.json').read_text())[label]
    actions[label] = {'revision': revision, 'kind': 'test', 'expected_exit': 0,
        'argv': ['python3', '-B', str(W / 'run_limb_benefit_slot.py'),
                 str(executables[mode].relative_to(R)), original['argv'][-1]],
        'limits': original['limits'], 'software_provider': True,
        'frozen_inputs': provider + resources + [str(executables[mode])] +
            [str(W / name) for name in ['software_contracts.py', 'compare_software.py', 'run_sequence.py',
             'run_limb_benefit_slot.py', 'baseline.json', 'build-and-arrays.json',
             'expected-readback-layout.json', 'retention-gate-frozen.json']] + [str(Q / 'audited-observer.cpp')],
        'scope': 'Corrected original audited observer ABBA; original full numerical workload,54 actual command joins and20+2 gates. First deficient rawcollection retained separately without verdict. No frame/full scientific/release qualification.'}
write(W / 'actions.json', actions)
freeze = json.loads((P / 'retention-gate-frozen.json').read_text())
freeze.update(candidate_revision=revision, baseline_revision='4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
    executables={mode: {'path': str(p.relative_to(R)), 'source_revision': revision if mode == 'candidate' else '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
        'bytes': p.stat().st_size, 'sha256': sha(p), 'observer_workspace_revision': revision,
        'observer_source_sha256': sha(Q / 'audited-observer.cpp')} for mode, p in executables.items()},
    actions_sha256=sha(W / 'actions.json'), scope='One admissible corrected original audited ABBA. Exact unchanged original reader/comparator/sequence/criteria, only executable allowlist adapted. Baseline4b product with original observer overlay; candidate8a product with same overlay. No initial deficient costverdict reused.')
assert all(freeze[key] == sha(W / name) for key, name in [('reader_sha256','software_contracts.py'),
    ('comparator_sha256','compare_software.py'),('launcher_sha256','run_sequence.py')])
write(W / 'retention-gate-frozen.json', freeze)
write(Q / 'audited-product-check.json', {'maps': maps, 'executables': freeze['executables'],
    'expected_readback_layout_sha256': sha(W / 'expected-readback-layout.json'),
    'original_reader_comparator_sequence_and20_plus2_criteria_exact': True})
shutil.copyfile(Q / 'actions.json', Q / 'actions-before-audited-witness.json')
parent_actions = json.loads((Q / 'actions.json').read_text())
parent_actions['audited-control-8'] = {'revision': revision, 'kind': 'test', 'expected_exit': 0,
    'argv': [str(executables['candidate']), '--gtest_filter=RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults'],
    'limits': actions['b1']['limits'], 'software_provider': True,
    'frozen_inputs': [str(executables['candidate'])] + provider + resources +
        [str(W / name) for name in ['software_contracts.py', 'build-and-arrays.json', 'expected-readback-layout.json']],
    'scope': 'Repaired original last correctness observer prerequisite; enforce32 complete digests/51 available markers/54 actual command joins and original inclusive API timing schema before correctedABBA. No performance verdict.'}
write(Q / 'actions.json', parent_actions)
print(json.dumps({'prepared': True, 'original_gates_unchanged': True, 'executables': freeze['executables']}))
