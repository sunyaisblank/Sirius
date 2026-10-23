"""Identify incremental outputs before the one owned compiler invocation."""
from pathlib import Path
import hashlib, json, re, shlex, subprocess

R = Path.cwd(); W = Path(__file__).resolve().parent
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
B = R / 'bin/linux-gcc/lavapipe-acqrel-ordering'
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes': p.stat().st_size,
    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
review = doc(W / 'review-admission.json'); assert review['passed']
assert pin(W / 'provider.patch') == review['prototype_patch']
assert pin(W / 'source-contract.json') == review['source_contract']
for key, value in review['source_files'].items(): assert pin(S / key) == value, key
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
dry = subprocess.check_output(['ninja', '-C', str(B), '-n', '-v',
    'src/gallium/targets/lavapipe/libvulkan_lvp.so'], text=True)
commands = []
for line in dry.splitlines():
    match = re.match(r'^\[\d+/\d+\] (.*)$', line)
    if match: commands.append(match[1])
assert commands, 'No new compiler work identified; stop before build.'
mutable = set()
for command in commands:
    words = shlex.split(command); outputs = set()
    for i, word in enumerate(words[:-1]):
        if word in {'-o', '--output', 'csr', 'csrD', 'csrDT', 'crs', 'crsD', 'rcs', 'rcsD'}:
            outputs.add(words[i + 1])
    assert outputs, ('Review unsupported incremental command before compilation', command)
    for output in outputs:
        p = (B / output).resolve()
        assert p.is_relative_to(B), p
        mutable.add(str(p.relative_to(R)))
target = B / 'src/gallium/targets/lavapipe/libvulkan_lvp.so'
assert str(target.relative_to(R)) in mutable
reused = {}
for p in B.rglob('*'):
    if not p.is_file(): continue
    key = str(p.relative_to(R))
    if key in mutable or p.name in {'.ninja_log', '.ninja_deps', '.ninja_lock'}: continue
    reused[key] = pin(p)
plan = {'source_independently_reviewed': True, 'source_review': pin(W / 'review-admission.json'),
    'previous_target': pin(target), 'dry_run_commands': commands,
    'mutable_outputs': sorted(mutable), 'reused_build_inputs': reused,
    'new_Meson_configuration': False, 'new_Sirius_CPP_or_shader_compilations': 0}
p = W / 'incremental-plan.json'; assert not p.exists()
p.write_text(json.dumps(plan, indent=2) + '\n')
(W / 'incremental-dry-run.txt').write_text(dry)
inputs = {str(p.relative_to(R)) for p in S.rglob('*') if p.is_file()}
inputs |= set(reused)
inputs |= {str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
inputs |= {'out/native-cpu-row-ownership/cold-frame/run_owned.py',
    'out/native-cpu-row-ownership/cold-frame/owned_guard.py', '/usr/bin/ninja', '/usr/bin/gcc-14',
    '/usr/bin/g++-14', '/usr/bin/ar', '/usr/bin/ld', '/usr/lib/llvm-21/bin/llvm-config',
    '/usr/lib/llvm-21/lib/libLLVM.so.21.1', str(R / 'bin/toolchains/python/bin/python3'),
    str(R / 'bin/toolchains/python/bin/meson')}
env = {'CC': '/usr/bin/gcc-14', 'CXX': '/usr/bin/g++-14',
    'PATH': str(R / 'bin/toolchains/python/bin') + ':/usr/lib/llvm-21/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin',
    'LLVM_CONFIG': '/usr/lib/llvm-21/bin/llvm-config', 'PYTHONDONTWRITEBYTECODE': '1'}
action = {'revision': revision, 'kind': 'build', 'expected_exit': 0,
    'limits': {'seconds': 180, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'argv': ['/usr/bin/python3', '-B', str(W / 'build_incremental.py')],
    'frozen_inputs': sorted(inputs), 'software_provider': False, 'build_environment': env,
    'scope': 'One incremental private provider build from independently reviewed additional exact E/DP arithmetic-module hashes with unchanged row-owner implementation. Existing original configuration and sealed reusable objects retained; only predeclared dirty outputs and Ninja logs/deps may change. Exact rejected acquire-release modification is absent. No global installation, new Sirius/probe/shader compilation, runtime/performance/science/release acceptance.'}
p = W / 'actions.json'; assert not p.exists()
p.write_text(json.dumps({'build': action}, indent=2) + '\n')
print(json.dumps({'source': revision, 'pending_commands': len(commands),
    'mutable_outputs': len(mutable), 'sealed_reused_build_inputs': len(reused),
    'original_config_reused': True}))
