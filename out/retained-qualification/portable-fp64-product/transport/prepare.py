#!/usr/bin/env python3
"""Prepare only: extract unchanged contract and link preset production code."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
STAGES = HERE.parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def run(args):
    r = subprocess.run(list(map(str,args)),cwd=ROOT,capture_output=True,text=True)
    if r.returncode: raise RuntimeError(f'{args}\n{r.stdout}\n{r.stderr}')
    return r.stdout

stage_report = json.loads((STAGES/'stage-artifact-report.json').read_text())
assert stage_report['pass'] and stage_report['dispatches'] == 0
(HERE/'modules').mkdir(exist_ok=True)
original_modules = {}
candidate_modules = {}
for stage in stage_report['stages']:
    stem = {'RayCamera':'ray_camera'}.get(stage['kind'],stage['kind'].lower())
    for mode in stage['modes']:
        source = STAGES/mode['path']
        assert sha(source) == mode['sha256']
        copied = HERE/'modules'/source.name
        shutil.copyfile(source,copied)
        candidate_modules[source.name] = sha(copied)
        suffix = '_portable_fp64' if mode['wide_product_definition'] else '_portable'
        original = ROOT/'bin/linux-gcc/src/sirius/backend/retained'/f'retained_{stem}{suffix}.spv'
        original_modules[f'{"fp64" if mode["wide_product_definition"] else "fp32"}:{stage["kind"]}'] = {
            'path':str(original.relative_to(ROOT)),'sha256':sha(original)}
shutil.copyfile(ROOT/'out/retained-qualification/portable-fp64-scalar/word_io.h', HERE/'word_io.h')
test_source = ROOT/'tests/backend/retained_compute_test.cpp'
text = test_source.read_text()
def block(marker):
    start=text.index(marker); opening=text.index('{',start); end=text.index('\n}\n',opening)+2
    return text[start:end], text[opening+1:end-1]
step_agrees,_=block('::testing::AssertionResult StepAgrees(')
_,body=block('TEST_F(RetainedComputeTest, JointRkStagesRetainCriticalIncrementsAndEmbeddedError)')
(HERE/'step_agrees_exact.txt').write_text(step_agrees)
(HERE/'transport_body_exact.txt').write_text(body)
template=(HERE/'runner_template.cpp').read_text()
runner=template.replace('// STEP_AGREES_EXACT_INSERT',step_agrees).replace('// TRANSPORT_BODY_EXACT_INSERT',body)
assert runner.count(step_agrees)==1 and runner.count(body)==1
(HERE/'runner.cpp').write_text(runner)
frozen_tests=HERE/'frozen-tests/support/retained_transport'
frozen_tests.mkdir(parents=True,exist_ok=True)
fixture_sources=[ROOT/'tests/support/retained_transport'/n for n in ['reference_cases.h','reference_cases.json','reference.py','README.md']]
for source in fixture_sources: shutil.copyfile(source,frozen_tests/source.name)
libs=[ROOT/'bin/linux-gcc/src/sirius/backend/libsirius_backend.a',ROOT/'bin/linux-gcc/src/sirius/core/libsirius_core.a',
      ROOT/'bin/linux-gcc/src/sirius/libsirius_base.a',ROOT/'bin/linux-gcc/lib/libgtest.a']
linked={str(p.relative_to(ROOT)):sha(p) for p in libs}
command=['/usr/bin/g++-14','-std=c++2c','-O3','-DNDEBUG','-Wall','-Wextra','-Wpedantic','-Werror',
         '-fno-fast-math','-ffp-contract=off','-DSIRIUS_CONTRACT_MODE=2','-DSIRIUS_HAS_RETAINED_COMPUTE=1',
         '-DSIRIUS_RETAINED_TESTS_AVAILABLE=1','-I',ROOT/'src','-I',HERE,'-I',HERE/'frozen-tests',
         '-I',ROOT/'bin/linux-gcc/src/sirius/backend/retained','-isystem',ROOT/'bin/linux-gcc/_deps/googletest-src/googletest/include',
         HERE/'runner.cpp','-o',HERE/'transport_runner',*libs,'-lvulkan','-ldl','-pthread']
run(command)
assert all(sha(ROOT/n)==h for n,h in linked.items()),'library changed during link'
preflight=run([HERE/'transport_runner',HERE/'modules',HERE/'evidence','--gtest_list_tests'])
assert preflight.count('JointRkStagesRetainCriticalIncrementsAndEmbeddedError/')==2
(HERE/'preflight.txt').write_text(preflight)
bound_production=[ROOT/'src/sirius/backend/retained_compute.h',ROOT/'src/sirius/backend/retained_compute.cpp',
                  ROOT/'src/sirius/backend/device.h',ROOT/'src/sirius/backend/vulkan/vulkan_device.cpp',
                  ROOT/'src/sirius/backend/vulkan/vulkan_device.h',ROOT/'src/sirius/core/twofold.h',test_source,
                  ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h',*fixture_sources]
report={'scope':'finite default-llvmpipe retained transport scalar/product-assist feasibility; not production/full-stage qualification',
        'source_head':run(['git','rev-parse','HEAD']).strip(),'tracked_working_tree_context':run(['git','status','--short','--untracked-files=no']),
        'stage_report_sha256':sha(STAGES/'stage-artifact-report.json'),
        'candidate_modules':candidate_modules,'original_modules':original_modules,'link_library_sha256':linked,
        'production_sources_sha256':{str(p.relative_to(ROOT)):sha(p) for p in bound_production},
        'command':list(map(str,command)),
        'contract':{'fixture_rows':15,'scientific_values_per_row':160,'expected_stages_per_row':7,
                    'group_relative_error_bound':'1e-11 unchanged',
                    'enclosure':'radius + independent precision gap + unchanged Twofold rounding bound',
                    'two_term_and_narrowed_negative_controls':28,'flat_embedded_error':'exact zero and equal fourth/fifth centers',
                    'capacity':24,'explicit_device_bytes_limit':8*1024*1024},
        'step_agrees_exact_sha256':sha(HERE/'step_agrees_exact.txt'),'transport_body_exact_sha256':sha(HERE/'transport_body_exact.txt'),
        'preflight':{'device_calls':0,'dispatches':0,'tests':preflight},
        'bound_diagnostic_files':{str(p.relative_to(HERE)):sha(p) for p in [HERE/'runner_template.cpp',HERE/'runner.cpp',
                     HERE/'word_io.h',HERE/'transport_runner',HERE/'prepare.py',*[HERE/'modules'/n for n in candidate_modules],
                     *list(frozen_tests.iterdir()),HERE/'step_agrees_exact.txt',HERE/'transport_body_exact.txt']}}
(HERE/'preparation-report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'prepared':True,'executable_sha256':sha(HERE/'transport_runner'),'device_calls':0,'dispatches':0,'tests':2}))
