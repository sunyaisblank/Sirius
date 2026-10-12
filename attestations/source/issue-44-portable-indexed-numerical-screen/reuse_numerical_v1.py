"""Require exact archived baseline/corpus/host joins before candidate execution."""
from pathlib import Path
import hashlib
import json
import resource
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024 * 1024, 64 * 1024 * 1024))
root = Path.cwd()
work = Path(__file__).resolve().parent
build = root / 'bin/linux-gcc/portable-indexed-limb-accumulation'
archive = root / 'attestations/software-vulkan/661e157/portable-limb-accumulation/files'
old = archive / 'out/portable-execution-attribution'
oldbuild = archive / 'bin/linux-gcc/portable-execution-attribution/limb-screen'
records = []
def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
def join(path, expected, reason):
    value = identity(path)
    assert value == {key: expected[key] for key in value}, (path, reason)
    records.append({'path': str(path.relative_to(root)), **value, 'reason': reason})
def receipt(case):
    folder = old / case
    owner = json.loads((folder / 'owner.json').read_text())
    before = json.loads((folder / 'inputs-before.json').read_text())
    after = json.loads((folder / 'inputs-after.json').read_text())
    assert owner['passed'] and owner['status'] == 'terminal' and not owner['cleanup_errors']
    assert not owner['remaining_owned_processes'] and owner['whole_inputs_unchanged'] and before == after
    return before
host = receipt('compile-limb-arithmetic-host')
for name in ('limb_arithmetic_gpu.cpp', 'triple_word_io.h'):
    join(old / name, host['out/portable-execution-attribution/' + name], 'original host compilation input')
    assert (work / name).read_bytes() == (old / name).read_bytes()
for name in ('src/sirius/backend/vulkan/vulkan_device.cpp', 'src/sirius/backend/vulkan/vulkan_device.h',
             'src/sirius/backend/device.h', 'src/sirius/base/error.h',
             'src/sirius/backend/vulkan/vulkan_portability.h', 'src/sirius/base/contracts.h',
             'bin/linux-gcc/src/sirius/libsirius_base.a', '/usr/lib/x86_64-linux-gnu/libvulkan.so'):
    join(root / name, host[name], 'unchanged compiled host/provider dependency')
for dep in ('limb-arithmetic-gpu.d', 'vulkan-device-v2.d'):
    text = (oldbuild / dep).read_text().replace('\\\n', '')
    for entry in text.split(':', 1)[1].split():
        path = Path(entry).resolve()
        key = str(path.relative_to(root))
        actual = archive / key if key.startswith('out/portable-execution-attribution/') else path
        join(actual, host[key], 'complete recorded non-system host/provider dependency')
recipe = json.loads((work / 'short-input-preset-recipe.json').read_text())
commands = json.loads((root / 'bin/linux-gcc/compile_commands.json').read_text())
assert recipe['recorded_CMake_adapter_command'] in commands
fixture = json.loads((old / 'limb-arithmetic-fixture.json').read_text())
assert fixture['records'] == 12763 and fixture['complete_words_each'] == 102104
assert fixture['requested_device_bytes'] <= 8 * 1024 * 1024
join(old / 'limb-arithmetic-input.bin', fixture['input'], 'unchanged corpus')
outputs = {}
for mode in ('integer', 'normal', 'normal_transport', 'normal_endpoint'):
    compiled = receipt('compile-limb-arithmetic-baseline-' + mode)
    for name in ('retained_scalar.slang', 'retained_pair.slang', 'retained_triple.slang', 'portable_binary32.h'):
        key = 'bin/linux-gcc/portable-execution-attribution/limb-screen/baseline-kernels/' + name
        join(root / 'src/sirius/kernels' / name, compiled[key], 'current baseline equals consumed source')
        join(oldbuild / 'baseline-kernels' / name, compiled[key], 'archived consumed baseline source')
    for key in compiled:
        if key.startswith(('/opt/slang/', '/usr/bin/spirv-opt')) or key.startswith('bin/linux-gcc/qualification-tools/') or key == 'scripts/build-retained-kernels.py':
            join(root / key, compiled[key], 'unchanged original compilation tool/recipe')
    consumed = receipt('gpu-limb-arithmetic-baseline-' + mode)
    for name in ('limb_arithmetic_gpu', 'limb-arithmetic-gpu.d', 'limb-arithmetic-gpu.map',
                 'vulkan-device-v2.o', 'vulkan-device-v2.d', 'arithmetic-baseline-' + mode + '.spv'):
        key = 'bin/linux-gcc/portable-execution-attribution/limb-screen/' + name
        join(oldbuild / name, consumed[key], 'baseline dispatch actual host/module input')
    for key in ('/usr/lib/x86_64-linux-gnu/libvulkan.so', '/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so',
                '/usr/lib/x86_64-linux-gnu/libLLVM.so.20.1', '/usr/share/vulkan/icd.d/lvp_icd.json',
                'tests/support/portable_binary32/reference.py'):
        join(root / key, consumed[key], 'unchanged provider/reference')
    verification = receipt('verify-limb-arithmetic-baseline-' + mode)
    output = old / ('gpu-limb-arithmetic-baseline-' + mode) / 'actual.bin'
    join(output, verification['out/portable-execution-attribution/gpu-limb-arithmetic-baseline-' + mode + '/actual.bin'], 'verified complete baseline output')
    report = json.loads((old / ('limb-arithmetic-check-baseline-' + mode + '.json')).read_text())
    assert report['passed'] and report['records'] == 12763
    assert identity(output) == {'bytes': report['output_bytes'], 'sha256': report['output_sha256']}
    assert report['all_valid_results_within_mathematical_preconditions_enclose_independent_exact_input_intervals']
    outputs[mode] = str(output.relative_to(root))
for name in ('limb-arithmetic-fixture.json', 'limb-arithmetic-input.bin'):
    assert not (work / name).exists()
    shutil.copy2(old / name, work / name)
for name in ('limb_arithmetic_gpu', 'limb-arithmetic-gpu.d', 'limb-arithmetic-gpu.map'):
    assert not (build / name).exists()
    shutil.copy2(oldbuild / name, build / name)
result = {'passed': True, 'exact_joins': records, 'baseline_outputs': outputs,
          'fresh_baseline_dispatches': 0, 'new_host_compilations': 0,
          'scope': 'Exact unchanged baseline source/corpus/provider/reference and consumed host/module/output reuse. Candidate arithmetic remains unexecuted; no current runtime benefit or qualification claim.'}
(work / 'numerical-reuse.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'passed': True, 'joins': len(records), 'baseline_modes': len(outputs)}))
