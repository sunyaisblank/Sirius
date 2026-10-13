"""Freeze the original eight control actions after linked-producer acceptance."""
from pathlib import Path
import hashlib, json, shutil

R = Path.cwd().resolve()
W = Path(__file__).resolve().parent
T = R / 'bin/linux-gcc/guarded-product-high-controls'
acceptance = json.loads((W / 'producer-acceptance.json').read_text())
assert acceptance['passed']
executable = T / 'original_candidate_controls'
assert hashlib.sha256(executable.read_bytes()).hexdigest() == acceptance['candidate_executable']['sha256']
prior = R / 'attestations/software-vulkan/41fe5f2/portable-denominator-preparation/files/out/portable-denominator-preparation/original-controls/actions.json'
original = json.loads(prior.read_text())
baseline_xml = R / 'attestations/software-vulkan/8a27321/vulkan-descriptor-reuse/files/out/vulkan-descriptor-reuse/audited-benefit/a1/gtest.xml'
shutil.copyfile(baseline_xml, W / 'original-baseline-timestamps.xml')
config = W / 'actions.json'
data = config.read_bytes()
versions = W / 'input-versions'
versions.mkdir(exist_ok=True)
(versions / hashlib.sha256(data).hexdigest()).write_bytes(data)
(W / 'actions-before-runtime.json').write_bytes(data)
actions = json.loads(data)
provider = ['/usr/share/vulkan/icd.d/lvp_icd.json', '/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so',
            '/usr/lib/x86_64-linux-gnu/libvulkan.so', '/usr/lib/x86_64-linux-gnu/libLLVM.so.20.1',
            '/lib/x86_64-linux-gnu/libc.so.6', '/lib/x86_64-linux-gnu/libm.so.6',
            '/lib/x86_64-linux-gnu/libstdc++.so.6', '/lib/x86_64-linux-gnu/libgcc_s.so.1',
            '/lib64/ld-linux-x86-64.so.2']
shared = list(actions['link-candidate']['frozen_inputs'])
shared += provider + [str(p) for p in [executable, T / 'original_baseline_controls', prior, baseline_xml,
           W / 'prepare_runtime.py', W / 'verify_producer.py', W / 'producer-acceptance.json',
           W / 'original-baseline-timestamps.xml', W / 'protected-before.json', W / 'input-seal-repair.json']]
for label in ['compile-compute-sealed', 'compile-observer', 'link-candidate']:
    shared += [str(W / label / name) for name in ['owner.json', 'inputs-before.json', 'inputs-after.json']]
shared += [str(p) for p in (W / 'contract').iterdir() if p.is_file()]
shared += [str(p) for p in (T / 'resources').rglob('*') if p.is_file()]
shared = list(dict.fromkeys(shared))
assert all(Path(p).is_file() for p in shared)
for index in range(1, 9):
    old = original['control-' + str(index)]
    label = 'control-' + str(index)
    assert label not in actions
    actions[label] = {'revision': actions['link-candidate']['revision'], 'kind': 'test', 'expected_exit': 0,
        'argv': [str(executable), old['argv'][1]], 'limits': old['limits'], 'software_provider': True,
        'frozen_inputs': shared,
        'scope': 'Unchanged original retained integrated control on matched private guarded-product candidate and pinned actual software Vulkan provider. Full original assertions, no tolerance/fixture/route weakening. No benefit or release verdict.'}
config.write_text(json.dumps(actions, indent=2) + '\n')
(W / 'runtime-preparation.json').write_text(json.dumps({'filters': [actions['control-'+str(i)]['argv'][1] for i in range(1,9)],
    'shared_input_count': len(shared), 'original_limits_unchanged': True,
    'producer_acceptance_sha256': hashlib.sha256((W / 'producer-acceptance.json').read_bytes()).hexdigest(),
    'baseline_timestamp_xml_sha256': hashlib.sha256(baseline_xml.read_bytes()).hexdigest(),
    'scope': 'Eight original integrated prerequisites; not a performance ABBA.'}, indent=2) + '\n')
print(json.dumps({'prepared': True, 'controls': 8, 'shared_inputs': len(shared), 'candidate': acceptance['candidate_executable']}))
