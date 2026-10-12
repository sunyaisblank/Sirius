from pathlib import Path
import json, shutil, subprocess

ROOT = Path.cwd()
WORK = ROOT / 'out/portable-stage-seven-owner-construction'
PHASE = WORK / 'original-controls-link'
BIN = ROOT / 'bin/linux-gcc/portable-stage-seven-owner'
PHASE.mkdir()
revision = subprocess.run(['git','rev-parse','HEAD'],capture_output=True,text=True,check=True).stdout.strip()
libraries = ['src/sirius/backend/libsirius_backend_cpu.a', 'lib/libgtest_main.a',
             'src/sirius/backend/libsirius_backend.a', 'src/sirius/render/libsirius_render.a',
             'lib/libgtest.a', 'src/sirius/backend/libsirius_backend_cpu.a',
             'src/sirius/backend/libsirius_backend.a', '/usr/lib/x86_64-linux-gnu/libvulkan.so',
             '-ldl', 'src/sirius/core/libsirius_core.a', 'src/sirius/libsirius_base.a']
libraries = [str(ROOT/'bin/linux-gcc'/f) if not f.startswith(('/', '-')) else f for f in libraries]
actions = {}
for producer, objects in {
    'baseline': [BIN/'original_baseline_test.o', ROOT/'bin/linux-gcc/tests/backend/CMakeFiles/sirius_backend_tests.dir/retained_dopri_test.cpp.o'],
    'candidate': [BIN/'original_candidate_test.o', BIN/'original_candidate_dopri.o', BIN/'compute.o', BIN/'integrator.o', BIN/'executor.o'],
}.items():
    executable = BIN / ('original_' + producer + '_controls')
    argv = ['/usr/bin/g++-14', '-O3', '-DNDEBUG', *map(str,objects), *libraries,
            '-Wl,-Map,' + str(PHASE/(producer+'.map')), '-o', str(executable)]
    frozen = [str(p.relative_to(ROOT)) for p in (WORK/'original-controls-preparation').rglob('*') if p.is_file() and not any(x in p.parts for x in ('compile-baseline-observer','compile-candidate-observer','compile-candidate-dopri','compile-candidate-executor'))]
    frozen += [str(p.relative_to(ROOT)) for p in objects]
    frozen += [f for f in libraries if not f.startswith('-')]
    frozen += ['bin/linux-gcc/build.ninja', 'bin/linux-gcc/compile_commands.json',
               'out/portable-stage-seven-owner-construction/prepare_original_links.py',
               '/usr/bin/g++-14', '/usr/bin/ld', '/usr/lib/x86_64-linux-gnu/libvulkan.so.1']
    actions['link-'+producer] = {'revision':revision,'kind':'diagnostic','expected_exit':0,
        'argv':argv,'limits':{'seconds':300,'rss_kib':4194304,'reserve_kib':2097152},
        'frozen_inputs':frozen,'scope': 'Diagnostic original-control link with configured preset libraries and a map; baseline9fbd source/ABI, candidate87 compute/integrator/executor; no full-product or qualification claim.'}
(PHASE/'actions.json').write_text(json.dumps(actions,indent=2)+'\n')
for file in ['run_owned.py','owned_guard.py']:
    shutil.copyfile(WORK/'runtime-stream-v2'/file, PHASE/file)
print(json.dumps({'phase':str(PHASE.relative_to(ROOT)),'actions':list(actions)},indent=2))
