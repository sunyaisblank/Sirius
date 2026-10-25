"""Seal one diagnostic build, then one component run; execute neither."""
from pathlib import Path
import hashlib
import json
import os
import re
import struct
import subprocess
import sys
from elf_proof import pin

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-accumulate-compatibility'
PLAN = R / 'out/native-cpu-accumulate-inline/coverage-plan.json'
RECEIPT = R / 'out/native-cpu-accumulate-inline/ispc-inline-compile-result.json'
LIMITS = {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2); stream.write('\n')


def owner_ok(path):
    d = json.loads(path.read_text())
    assert d['status'] == 'terminal' and d['passed'] and d['source_unchanged'] and d['whole_inputs_unchanged']
    assert not d['remaining_owned_processes'] and not d['cleanup_errors']
    assert d['observed_births'] and all(v['absent'] for v in d['observed_births'])


mode, = sys.argv[1:]
assert mode in ('compile', 'run')
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not any(os.environ.get(k) for k in ('CPATH', 'CPLUS_INCLUDE_PATH', 'LIBRARY_PATH', 'GCC_EXEC_PREFIX', 'COMPILER_PATH'))
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
plan = json.loads(PLAN.read_text()); compiled = json.loads(RECEIPT.read_text())
assert compiled['all_objects_compiled'] and compiled['all_prepared_source_bytes_unchanged']
assert compiled['new_compilations'] == 3 and compiled['reused_compiled_objects'] == 0
assert [v['stage'] for v in compiled['objects']] == ['transport', 'endpoint', 'dopri_phase']
new = [v['object'] for v in compiled['objects']]
for stage, row, record in zip(('transport', 'endpoint', 'dopri_phase'), new, compiled['objects']):
    assert row['path'] == f'bin/linux-gcc/native-cpu-accumulate-inline/objects/{stage}.o'
    assert record['exit'] == 0 and record['argv_executed'] and record['status'] == 'terminal'
    assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}
q_path = R / 'out/native-cpu-transport-worker-occupancy/provider/source-contract.json'
q = json.loads(q_path.read_text())['native_objects']; assert len(q) == 8
objects = q[:4] + new + q[7:]; assert len(objects) == 8
baseline = plan['standalone_link']['baseline_object']
assert pin(R / baseline['path']) == {k: baseline[k] for k in ('bytes', 'sha256')}
inputs = {**compiled['source_inputs'], **plan['source_inputs']}
for item in objects:
    assert pin(R / item['path']) == {k: item[k] for k in ('bytes', 'sha256')}
    inputs[item['path']] = {k: item[k] for k in ('bytes', 'sha256')}
inputs[baseline['path']] = {k: baseline[k] for k in ('bytes', 'sha256')}
for p in (PLAN, RECEIPT, q_path, *W.iterdir()):
    if p.is_file() and p.name not in ('prepared.json', 'actions.json'):
        inputs[str(p.relative_to(R))] = pin(p)
if mode == 'compile':
    assert not (W / 'prepared.json').exists() and not (W / 'actions.json').exists() and not B.exists()
    adapter_path = R / 'out/native-cpu-simd-packet-policy/adapter-compile-result.json'
    adapter = json.loads(adapter_path.read_text())
    assert adapter['all_prepared_source_bytes_unchanged']
    item, = adapter['objects']; assert item['exit'] == 0
    assert item['object'] == objects[-1]
    inputs.update(adapter['source_inputs']); inputs[str(adapter_path.relative_to(R))] = pin(adapter_path)
    # Same GNU14's completed dependency evidence; no new preprocessor/probe run.
    prior = R / 'out/native-cpu-host-representation/probe'
    fp_compile = json.loads((prior / 'compile-result.json').read_text())
    fp_prepare = json.loads((prior / 'preparation.json').read_text())
    assert fp_compile['passed'] and fp_compile['exit'] == 0
    compiler_path = str(Path(item['argv'][0]).resolve(strict=True))
    assert pin(Path(compiler_path)) == fp_prepare['tool_dependencies'][compiler_path]
    assert any(Path(p).name == 'cc1plus' for p in fp_prepare['tool_dependencies'])
    assert any(Path(p).name == 'collect2' for p in fp_prepare['tool_dependencies'])
    inputs.update(fp_prepare['tool_dependencies'])
    inputs.update({p: value for p, value in fp_compile['headers'].items() if not p.endswith('.cpp')})
    for p in (prior / 'compile-result.json', prior / 'preparation.json'):
        inputs[str(p.relative_to(R))] = pin(p)
    flags = list(item['argv'])
    assert all(v in flags for v in ('-O2', '-march=x86-64', '-mfpmath=sse', '-fno-fast-math',
                                   '-fno-associative-math', '-fno-unsafe-math-optimizations',
                                   '-frounding-math', '-ffp-contract=off'))
    original_flags = list(flags)
    flags[flags.index('-c') + 1] = str(W / 'driver.cpp')
    flags[flags.index('-o') + 1] = str(B / 'driver.o')
    flags += ['-MD', '-MF', str(W / 'compile/driver.d')]
    link = [v for i, v in enumerate(original_flags) if v not in ('-c', '-o') and
            i not in (original_flags.index('-c') + 1, original_flags.index('-o') + 1)]
    link += [str(B / 'driver.o'), *[str(R / v['path']) for v in objects], str(B / 'baseline.o'),
             '-o', str(B / 'compatibility'), '-Wl,--no-relax', '-lm']
    mapping = dict(line.split() for line in (W / 'baseline-symbols.map').read_text().splitlines() if line)
    assert mapping == plan['standalone_link']['rename_only_global_defined_symbols']
    # Compare the actual cached canonical array to the small original wire array.
    canonical = (R / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h').read_text()
    body = re.search(r'kTransportProgram\s*\{\{(.*?)\}\};', canonical, re.S)
    assert body
    wire = (R / plan['inputs']['program']['wire_source']).read_text()
    wire_body = re.search(r'transport_program\[\]\s*=\s*\{(.*?)\};', wire, re.S)
    assert wire_body
    words = lambda text: [int(v, 0) for v in re.findall(r'\b(0x[0-9a-fA-F]+|[0-9]+)[uU]?\b', text)]
    a, b = words(body[1]), words(wire_body[1])
    assert a == b and len(a) == 15967
    assert hashlib.sha256(struct.pack('<15967I', *a)).hexdigest() == plan['inputs']['program']['sha256']
    for name in (flags[0], '/usr/bin/nm', '/usr/bin/objcopy', '/usr/bin/ld', '/usr/bin/as',
                 '/usr/bin/ldd', '/usr/bin/python3'):
        p = Path(name).resolve(strict=True); inputs[str(p)] = pin(p)
    for name, expected in inputs.items():
        assert pin(R / name) == {k: expected[k] for k in ('bytes', 'sha256')}, name
    B.mkdir(parents=True, exist_ok=False)
    prepared = {'revision': revision, 'source_inputs': inputs, 'plan': pin(PLAN), 'R_compile_receipt': pin(RECEIPT),
                'objects': objects, 'candidate_transport_object': new[0], 'original_baseline_object': baseline,
                'rename_mapping': mapping, 'rename_argv': ['/usr/bin/objcopy',
                    '--redefine-syms=' + str(W / 'baseline-symbols.map'), str(R / baseline['path']), str(B / 'baseline.o')],
                'compile_argv': flags, 'link_argv': link, 'original_I_compile_argv': original_flags,
                'link_scope': 'Same baseline-ISA/strict-FP GNU14 flags and libc/libm/C++ runtime. GNU --no-relax preserves opcodes for declared4B-relocation association; no math/compiler option changed.',
                'header_tool_provenance': 'Pre-frozen actual same-GNU14 cc1plus/collect2 and completed standard/intrinsic-header dependency pins from prior FP compile, plus exact fixture/wire/rows/native-bridge plan inputs. Reuse of dependency evidence only; no FP-matrix execution.',
                'baseline': str((B / 'baseline.o').relative_to(R)), 'driver': str((B / 'driver.o').relative_to(R)),
                'executable': str((B / 'compatibility').relative_to(R)), 'limits': LIMITS}
    write(W / 'prepared.json', prepared); actions = {}
else:
    owner_ok(W / 'compile/owner.json')
    build = json.loads((W / 'compile-result.json').read_text()); assert build['passed']
    prepared = json.loads((W / 'prepared.json').read_text())
    assert build['prepared'] == pin(W / 'prepared.json') and build['candidate_transport_object'] == new[0]
    inputs.update(prepared['source_inputs'])
    inputs.update(build['runtime_libraries'])
    for key in ('renamed_baseline_object', 'driver_object', 'executable', 'link_proof', 'rename_proof'):
        row = build[key]; inputs[row['path']] = {k: row[k] for k in ('bytes', 'sha256')}
    for p in (W / 'compile-result.json', W / 'compile/owner.json'):
        inputs[str(p.relative_to(R))] = pin(p)
    actions = json.loads((W / 'actions.json').read_text()); assert 'run' not in actions and not (W / 'run').exists()
inputs[str((W / 'prepared.json').relative_to(R))] = pin(W / 'prepared.json')
actions[mode] = {'revision': revision, 'kind': 'diagnostic', 'expected_exit': 0, 'software_provider': False,
                 'limits': LIMITS, 'frozen_inputs': sorted(inputs),
                 'argv': ['/usr/bin/python3', '-B', str(W / 'task.py'), mode, '{folder}'],
                 'scope': 'One separately owned R raw Transport compatibility ' + mode +
                 '; ALL38 original fixture inputs/187960 output words. One attempt, no timing/science/original-gate substitution.'}
if mode == 'compile':
    write(W / 'actions.json', actions)
else:
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'prepared': mode, 'executed': False, 'limits': LIMITS}))
