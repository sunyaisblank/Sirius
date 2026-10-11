"""Build two isolated current factories; no Vulkan workload or production write."""
from pathlib import Path
import hashlib, json, os, re, shlex, shutil, subprocess, sys, time
import payload_support as payload
sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[1]
BUILD = ROOT / 'bin/linux-gcc'
REV = 'b81e061fa081cdaed37eca3add950bd254bf296d'
CHANGED = {'kTransportPortableNormalSumShader'}
CASES = ['RetainedComputeTest.JointRkStagesRetainCriticalIncrementsAndEmbeddedError',
         'RetainedComputeTest.MixedIndependentTransportLayersPreserveWordsAndRefusal',
         'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults',
         'RetainedComputeTest.WideDeviceTimestampsPreserveOriginalIntervalResults']

def seal(path):
    path = Path(path).resolve(strict=True)
    return {'path': str(path), 'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}

def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')

def source():
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip() == REV
    assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT)
    assert subprocess.check_output(['git', 'rev-list', '--left-right', '--count', 'HEAD...@{upstream}'], cwd=ROOT).split() == [b'0', b'0']
    assert subprocess.check_output(['git', 'worktree', 'list', '--porcelain'], cwd=ROOT).count(b'worktree ') == 1

def main():
    source()
    compiled = json.loads((WORK / 'compile-result.json').read_text())
    assert compiled['pass'] and seal(WORK / 'literal-transport.spv') == compiled['module']
    assert json.loads((WORK / 'compiled-inspection.json').read_text())['pass']
    header = BUILD / 'src/sirius/backend/retained/retained_kernels.h'
    original_text = header.read_text()
    original = payload.arrays(original_text)
    candidate_text = payload.replace_payload(original_text, next(iter(CHANGED)), (WORK / 'literal-transport.spv').read_bytes())
    candidate = payload.arrays(candidate_text)
    assert {n for n in original if original[n] != candidate[n]} == CHANGED
    assert payload.ARRAY.sub('<ARRAY>', original_text) == payload.ARRAY.sub('<ARRAY>', candidate_text)
    canonical = (ROOT / 'src/sirius/backend/retained_compute.cpp').read_text()
    modules = WORK / 'modules'; modules.mkdir()
    for mode, text in [('baseline', original_text), ('candidate', candidate_text)]:
        namespace = payload.NAMESPACE + mode
        private = payload.replace_once(text, 'namespace sirius::backend::retained_program {', 'namespace sirius::backend::' + namespace + ' {')
        private = payload.replace_once(private, '// namespace sirius::backend::retained_program', '// namespace sirius::backend::' + namespace)
        (modules / (mode + '_retained_kernels.h')).write_text(private)
        factory = payload.replace_once(canonical, '#include "retained_kernels.h"', '#include "modules/' + mode + '_retained_kernels.h"')
        factory = payload.replace_once(factory, 'using namespace retained_program;', 'using namespace ' + namespace + ';')
        (WORK / (mode + '_retained_compute.cpp')).write_text(factory)
        assert factory.replace('#include "modules/' + mode + '_retained_kernels.h"', '#include "retained_kernels.h"').replace('using namespace ' + namespace + ';', 'using namespace retained_program;') == canonical
    entries = json.loads((BUILD / 'compile_commands.json').read_text())
    entry = next(e for e in entries if e['file'] == str(ROOT / 'src/sirius/backend/retained_compute.cpp'))
    base = shlex.split(entry['command']); compiler = Path(base[0]).resolve(strict=True)
    assert all(flag in base for flag in ['-Werror', '-std=c++2c', '-DSIRIUS_CONTRACT_MODE=2'])
    target = next(line for line in (BUILD / 'build.ninja').read_text().splitlines() if line.startswith('build ') and ': CXX_EXECUTABLE_LINKER__sirius_backend_tests_' in line)
    items = shlex.split(target.split(': ', 1)[1])[1:]
    objects = [BUILD / item for item in items[:items.index('|')]]
    assert len(objects) == 20 and all(p.suffix == '.o' for p in objects)
    libs = [BUILD / n for n in ['src/sirius/backend/libsirius_backend_cpu.a', 'src/sirius/backend/libsirius_backend.a', 'src/sirius/render/libsirius_render.a', 'src/sirius/core/libsirius_core.a', 'src/sirius/libsirius_base.a', 'lib/libgtest_main.a', 'lib/libgtest.a']]
    normal_names = subprocess.check_output(['nm', '-C', '--defined-only', next(p for p in objects if p.name == 'retained_compute_test.cpp.o')], text=True)
    normal = {name: data for name, data in original.items() if ' sirius::backend::retained_program::' + name + '\n' in normal_names}
    assert len(normal) == 36
    inputs = [ROOT / os.fsdecode(p) for p in subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).split(b'\0') if p]
    inputs += [*objects, *libs, compiler, header, BUILD / 'CMakeCache.txt', BUILD / 'compile_commands.json', BUILD / 'build.ninja', Path('/usr/lib/x86_64-linux-gnu/libvulkan.so').resolve(strict=True)]
    inputs += list((BUILD / 'src/sirius/backend/retained').glob('*.spv'))
    assert len(list((BUILD / 'src/sirius/backend/retained').glob('*.spv'))) == 33
    inputs += [WORK / name for name in ['prepare_hosts.py', 'payload_support.py', 'compile-before.json', 'compile-result.json', 'compiled-inspection.json', 'literal-transport.spv']]
    inputs += list(modules.iterdir()) + [WORK / (mode + '_retained_compute.cpp') for mode in ('baseline', 'candidate')]
    for option in ['-print-prog-name=cc1plus', '-print-prog-name=collect2', '-print-prog-name=ld']:
        name = subprocess.check_output([compiler, option], text=True).strip()
        inputs.append(Path(name if '/' in name else shutil.which(name)).resolve(strict=True))
    commands = []
    def run(label, argv, cwd=ROOT):
        argv = [str(x) for x in argv]; started = time.monotonic()
        with (WORK / (label + '.stdout')).open('wb') as out, (WORK / (label + '.stderr')).open('wb') as err:
            result = subprocess.run(argv, cwd=cwd, stdout=out, stderr=err, check=False)
        commands.append({'label': label, 'argv': argv, 'cwd': str(cwd), 'exit': result.returncode, 'seconds': time.monotonic() - started})
        write(WORK / 'host-commands.json', commands)
        assert result.returncode == 0, label
    plans = {}
    for mode in ('baseline', 'candidate'):
        obj = WORK / (mode + '_retained_compute.o'); assert not obj.exists()
        argv = base.copy(); argv[0] = str(compiler); argv[argv.index('-o') + 1] = str(obj); argv[-1] = str(WORK / (mode + '_retained_compute.cpp'))
        plans[mode] = argv
        dependency = WORK / (mode + '.d'); scan = argv.copy(); scan[scan.index('-o') + 1] = str(dependency)
        scan[-1:-1] = ['-M', '-MF', str(dependency)]
        run(mode + '-dependencies', scan, Path(entry['directory']))
        names = shlex.split(dependency.read_text().replace('\\\n', ' ').split(':', 1)[1]); assert names
        inputs += [Path(n) if Path(n).is_absolute() else Path(entry['directory']) / n for n in names]
        inputs.append(dependency)
    before = {str(p.resolve(strict=True)): seal(p) for p in inputs}
    write(WORK / 'host-inputs-before.json', before)
    result = {'pass': False, 'device_workloads': 0, 'revision': REV, 'changed_arrays': sorted(CHANGED), 'unchanged_other_arrays': 39, 'factory_source_two_substitutions_only': True, 'artifacts': {}}
    try:
        for mode, arrays in [('baseline', original), ('candidate', candidate)]:
            run(mode + '-factory', plans[mode], Path(entry['directory']))
            exe = WORK / (mode + '-backend-tests'); assert not exe.exists()
            run(mode + '-link', [compiler, '-O3', '-DNDEBUG', WORK / (mode + '_retained_compute.o'), *objects, '-Wl,--start-group', *libs, '-Wl,--end-group', '-lvulkan', '-ldl', '-pthread', '-o', exe])
            private_spans = payload.elf_arrays(exe, arrays, payload.NAMESPACE + mode, normal)
            normal_spans = payload.elf_arrays(exe, normal, 'retained_program')
            assert len(private_spans) == 40 and len(normal_spans) == 36
            run(mode + '-inventory', [exe, '--gtest_list_tests', '--gtest_filter=' + ':'.join(CASES)])
            listed = (WORK / (mode + '-inventory.stdout')).read_text().splitlines()
            assert [line.strip() for line in listed if line.startswith('  ')] == [case.split('.')[1] for case in CASES]
            result['artifacts'][mode] = {'executable': seal(exe), 'factory': seal(WORK / (mode + '_retained_compute.o')), 'all40_private_spans': private_spans, 'original36_test_reference_spans': normal_spans}
        result['pass'] = True
    finally:
        after = {name: seal(name) for name in before}
        write(WORK / 'host-inputs-after.json', after)
        source(); result['all_inputs_unchanged'] = before == after
        result['input_count'] = len(before)
        result['pass'] = result['pass'] and result['all_inputs_unchanged']
        write(WORK / 'host-preparation.json', result)
    assert result['pass']
    print(json.dumps({'pass': result['pass'], 'inputs': len(before), 'private_arrays_each': 40, 'original_test_spans_each': 36, 'device_workloads': 0}))

if __name__ == '__main__':
    main()
