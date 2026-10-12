"""Freeze the unchanged original comparison for one mixed-composition trial."""
from pathlib import Path
import hashlib, json, resource, shutil, subprocess, sys
sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64*1024*1024, 64*1024*1024))
r = Path.cwd()
w = Path(__file__).resolve().parent
p = w / 'benefit'
b = r / 'bin/linux-gcc/portable-transport-limb-accumulation'
saved = r / 'attestations/software-vulkan/661e157/portable-limb-accumulation/files/out/portable-execution-attribution/benefit'
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert revision == 'd910f3c3b186bb5ed28a8102faf84b24698d679d'
assert subprocess.check_output(['git', 'status', '--porcelain'], text=True) == ''
for name in ['consumer-check.json', 'transport-generated-products-check.json',
             'independent-consumer-review.json', 'independent-integrated-control-review.json',
             'integrated-control-check.json']:
    assert json.loads((w/name).read_text())['passed'], name
p.mkdir()
for name in ['software_contracts.py', 'compare_software.py', 'run_sequence.py',
             'expected-readback-layout.json']:
    shutil.copy2(saved/name, p/name)
    assert sha(saved/name) == sha(p/name)
for name in ['run_owned.py', 'owned_guard.py', 'run_limb_benefit_slot.py']:
    shutil.copy2(w/name, p/name)
for source, destination in [('transport-baseline-arrays.json', 'baseline.json'),
                            ('transport-candidate-arrays.json', 'build-and-arrays.json')]:
    shutil.copy2(w/source, p/destination)
(p/'edits.json').write_text(json.dumps({
    'previous_corrected_reader_and_comparator_copied_byte_exact': True,
    'all20_numeric_and2_inclusive_API_conditions_unchanged': True,
    'previous_reader_sha256': sha(saved/'software_contracts.py'),
    'previous_comparator_sha256': sha(saved/'compare_software.py'),
    'new_reader_or_comparator_edits': []}, indent=2)+'\n')
executables = {}
for mode, source in [('baseline', '9fbd0e8b0d549a13476d066f2877dcd2cf9a94d9'), ('candidate', revision)]:
    exe = b / ('original_'+mode+'_controls')
    executables[mode] = {'path': str(exe.relative_to(r)), 'source_revision': source,
                         'bytes': exe.stat().st_size, 'sha256': sha(exe)}
frozen = [str(x.relative_to(r)) for x in p.iterdir() if x.is_file()]
frozen += [x['path'] for x in executables.values()]
frozen += [str(b/'products/retained_kernels.h'), str(r/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'),
           '/usr/lib/x86_64-linux-gnu/libvulkan.so', '/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so',
           '/usr/lib/x86_64-linux-gnu/libLLVM.so.20.1', '/usr/share/vulkan/icd.d/lvp_icd.json', str(__file__)]
frozen += [str(w/name) for name in ['consumer-check.json', 'transport-generated-products-check.json',
           'independent-consumer-review.json', 'finite-arithmetic-evidence-reuse.json',
           'independent-integrated-control-review.json', 'baseline-reuse.json', 'integrated-control-check.json',
           'independent-transport-source-review.json', 'independent-transport-products-review.json',
           'independent-consumer-control-owner-closure-review.json']]
frozen += [str(b/(name+suffix)) for name in ['candidate-observer', 'candidate-dopri', 'compute', 'integrator', 'trace']
           for suffix in ['.o', '.o.d']]
frozen += [str(b/name) for name in ['candidate.map', 'baseline.map', 'vulkan.o']]
assert all((r/name).is_file() for name in frozen)
actions = {}
for label in ['a1', 'b1', 'b2', 'a2']:
    mode = 'baseline' if label.startswith('a') else 'candidate'
    actions[label] = {'revision': revision, 'kind': 'test', 'expected_exit': 0,
        'argv': ['python3', '-B', str(p/'run_limb_benefit_slot.py'), executables[mode]['path'],
                 '--gtest_filter=RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults'],
        'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
        'software_provider': True, 'frozen_inputs': frozen+[str(p/'retention-gate-frozen.json')],
        'scope': 'One original default timestamp slot. Require all original 32 complete buffers, 51 ordered markers, 54 actual command joins, 20 performance and both inclusive API-cost conditions. No cost subtraction, retry or pooled order.'}
(p/'actions.json').write_text(json.dumps(actions, indent=2)+'\n')
freeze = {'candidate_revision': revision, 'baseline_revision': executables['baseline']['source_revision'],
    'executables': executables, 'sequence': ['a1', 'b1', 'b2', 'a2'], 'within_slot_repeats': 3,
    'complete_call_minimum_improvement_percent': 5, 'every_inclusive_stage_minimum_improvement_percent': -5,
    'all20_gates_required': True, 'observed_instance_api_minimum_improvement_percent': 0,
    'reader_sha256': sha(p/'software_contracts.py'), 'comparator_sha256': sha(p/'compare_software.py'),
    'launcher_sha256': sha(p/'run_sequence.py'), 'actions_sha256': sha(p/'actions.json'),
    'scope': 'Original comparison copied byte exact. One frozen ABBA, all original words, streams and costs required. Mixed-composition integration and consumer prerequisites pass; full qualification remains separate.'}
(p/'retention-gate-frozen.json').write_text(json.dumps(freeze, indent=2)+'\n')
print(json.dumps(freeze))
