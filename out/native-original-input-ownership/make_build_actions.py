from pathlib import Path
import json, re, shlex, subprocess, sys
R = Path.cwd(); W = Path(__file__).resolve().parent; C = R / 'bin/linux-gcc'
stage = sys.argv[1]; assert stage in ['modules-completion', 'matched']
inputs = {str(p.relative_to(R)) for p in W.rglob('*') if p.is_file()}
inputs |= {str(p.relative_to(R)) for p in (C / 'src/sirius/backend/retained').iterdir() if p.is_file()}
inputs.add('bin/linux-gcc/native-retained-ownership/modules/owned-dopri_phase.spv')
inputs.add('bin/linux-gcc/native-original-input-ownership/modules/owned-transport.spv')
inputs |= {'/usr/bin/python3', '/etc/ld.so.cache', 'bin/linux-gcc/CMakeCache.txt',
           'bin/linux-gcc/compile_commands.json', 'bin/linux-gcc/build.ninja'}
tools = ['/usr/bin/g++-14', '/usr/bin/ar', '/usr/bin/ninja', '/opt/slang/bin/slangc']
for line in (C / 'CMakeCache.txt').read_text().splitlines():
    if line.startswith('SIRIUS_SPIRV_') and ':FILEPATH=' in line:
        tools.append(line.split('=', 1)[1])
if stage == 'matched':
    assert json.loads((W / 'modules-completion/owner.json').read_text())['passed']
    lines = subprocess.check_output(['ninja', '-C', str(C), '-t', 'commands',
        'tests/backend/sirius_backend_tests'], text=True).splitlines()
    line, = [x for x in lines if ' -o tests/backend/sirius_backend_tests ' in x]
    for arg in shlex.split(line):
        p = Path(arg); p = p if p.is_absolute() else C / p
        if p.is_file(): inputs.add(str(p.relative_to(R)) if p.is_relative_to(R) else str(p))
    inputs.add('bin/linux-gcc/native-original-input-ownership/matched/audited-observer.o')
for name in tools:
    p = Path(name); inputs.add(str(p.resolve()))
    answer = subprocess.run(['ldd', str(p)], text=True, capture_output=True)
    inputs |= {str(Path(x).resolve()) for x in re.findall(r'(/[^\s()]+)', answer.stdout) if Path(x).is_file()}
a = {'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
     'kind': 'diagnostic', 'expected_exit': 0,
     'limits': {'seconds': 180, 'rss_kib': 4194304, 'reserve_kib': 2097152},
     'argv': ['/usr/bin/python3', '-B', str(W / ('build_modules.py' if stage == 'modules-completion' else 'build_matched.py'))],
     'frozen_inputs': sorted(inputs), 'software_provider': False,
     'scope': ('Reuse the exact already-compiled/validated Transport after a pure layout-reader correction; compile exactly one native32 Endpoint module and reuse the prior byte-exact DopriPhase module with the original builder, controls, programs, LocalSize64 and scalar algorithms; preserve all35 original modules. No Vulkan execution or qualification.' if stage == 'modules-completion' else 'Compile exactly one private three-module selector/object and relink the original audited observer with all other existing strict GCC objects; preserve every other route and archive member. No Vulkan execution or qualification.')}
p = W / 'actions.json'; actions = json.loads(p.read_text()) if p.exists() else {}
assert stage not in actions
actions[stage] = a; p.write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'stage': stage, 'inputs': len(inputs)}))
