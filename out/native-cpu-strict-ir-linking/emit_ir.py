"""Expose the existing DP compiler boundary before selecting any LTO tooling."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import re
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-strict-ir-linking'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


phase, = sys.argv[1:]
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    assert not B.exists() and not (W / 'actions.json').exists()
    parent = R / 'out/native-cpu-inline-private-bank/derivation.json'
    proof = json.loads(parent.read_text())
    command, = [row for row in proof['raw_compilations'] if row['stage'] == 'dopri_phase']
    argv = list(command['argv'])
    assert argv[-3:-1] == ['--emit-obj', '-o']
    argv[-3], argv[-1] = '--emit-llvm-text', str(B / 'dopri_phase.ll')
    inputs = dict(proof['source_inputs']) | dict(proof['derived'])
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    for path in [parent, *W.glob('*.py'),
                 R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
                 R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py']:
        inputs[str(path.relative_to(R))] = pin(path)
    write(W / 'plan.json', dict(argv=argv, source_inputs=inputs,
        parent_derivation=dict(path=str(parent.relative_to(R)), **pin(parent)),
        only_compiler_output_kind_and_path_changed=True))
    write(W / 'actions.json', {'emit': dict(
        revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=['/usr/bin/python3', '-B', str(W / 'emit_ir.py'), 'emit'],
        frozen_inputs=sorted(inputs) + [str((W / 'plan.json').relative_to(R))],
        scope='One unchanged smallest AS raw source emitted as LLVM text using the same pinned ISPC/target/O2/PIC/disable-fma/system-math flags. No new runtime, fixture, frame, benchmark, FP flag relaxation or tool acquisition. Actual FP vocabulary/attributes decide whether a separately justified constrained conversion is feasible.')})
    print(json.dumps(dict(prepared='emit', executed=False, frozen_inputs=len(inputs))))
elif phase == 'emit':
    plan = json.loads((W / 'plan.json').read_text())
    for name, expected in plan['source_inputs'].items():
        assert pin(R / name) == expected, name
    B.mkdir(exist_ok=False)
    subprocess.run(plan['argv'], check=True)
    output = B / 'dopri_phase.ll'
    text = output.read_text()
    operations = Counter()
    calls = Counter()
    flags = Counter()
    examples = {}
    fp_ops = {'fadd', 'fsub', 'fmul', 'fdiv', 'frem', 'fneg', 'fcmp', 'fpext',
              'fptrunc', 'fptosi', 'fptoui', 'sitofp', 'uitofp'}
    for line in text.splitlines():
        instruction = re.match(r'^\s+(?:%[^=]+ = )?([a-z]+)\b(.*)', line)
        if not instruction:
            continue
        op = instruction[1]
        if op in fp_ops:
            operations[op] += 1
            examples.setdefault(op, line.strip())
            flags.update(re.findall(r'\b(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b', line))
        if re.search(r'\bcall\b', line) and re.search(r'\b(?:float|double|half|fp128)\b', line):
            target = re.search(r'@([^ (]+)\(', line)
            if target:
                calls[target[1]] += 1
    attributes = dict(re.findall(r'^attributes #(\d+) = (.*)$', text, re.M))
    definitions = [line for line in text.splitlines() if line.startswith('define ')]
    strict = []
    for line in definitions:
        group = re.search(r'#(\d+)', line)
        if 'strictfp' in line or (group and 'strictfp' in attributes[group[1]]):
            strict.append(line)
    result = dict(emitted=True, output=dict(path=str(output.relative_to(R)), **pin(output)),
        argv=plan['argv'], definitions=len(definitions), strictfp_definitions=strict,
        ordinary_fp_operations=dict(operations), fp_instruction_flags=dict(flags),
        fp_examples=examples, floating_call_targets=dict(calls), attribute_groups=attributes,
        triple=re.search(r'^target triple = (.*)$', text, re.M)[1],
        data_layout=re.search(r'^target datalayout = (.*)$', text, re.M)[1],
        constrained_intrinsic_mentions=text.count('@llvm.experimental.constrained.'),
        scope='Compiler IR inventory only; no conversion/link/inlining, numerical equivalence, speed or adoption established. Attribute-only strictfp remains unsupported. Every FP operation/call requires a justified exact mapping before selecting a new pipeline.')
    write(W / 'inventory.json', result)
    print(json.dumps({key: result[key] for key in ('definitions', 'ordinary_fp_operations', 'fp_instruction_flags', 'floating_call_targets', 'constrained_intrinsic_mentions')}))
else:
    raise ValueError(phase)
