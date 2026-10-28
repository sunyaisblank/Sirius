"""Lower the structurally admitted DP copy once and inspect actual machine code."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-strict-ir-linking'
D = B / 'toolchain/bin'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,
                    sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


phase, = sys.argv[1:]
features = phase.endswith('-features')
base_phase = phase.removesuffix('-features') if features else phase
plan_path = W / ('machine-features-plan.json' if features else 'machine-plan.json')
case = 'machine-dp-features' if features else 'machine-dp'
target = B / ('dp-value.features.o' if features else 'dp-value.o')
if base_phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    owner = json.loads((W / 'value-optimize-joins/owner.json').read_text())
    assert owner['passed'] and owner['status'] == 'terminal'
    product = json.loads((W / 'value-optimize-joins-result.json').read_text())
    assert product['adapter_calls'] == product['packed_payload_allocas'] == product['packed_tag_allocas'] == 0
    assert product['core_definitions'] == 8
    paths = [*W.glob('*.py'), W / 'value_abi.cpp', W / 'value_abi_contract.json',
             W / 'scope_contract.json', W / 'value-optimize-joins-result.json',
             W / 'toolchain.json', R / product['output']['path'],
             R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
             R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py',
             Path('/usr/bin/objdump').resolve(), Path('/usr/bin/nm').resolve()]
    row = json.loads((W / 'toolchain.json').read_text())['tools']['llc']
    assert pin(R / row['resolved']) == {k: row[k] for k in ('bytes', 'sha256')}
    paths.append(R / row['resolved'])
    if features:
        failed = json.loads((W / 'machine-dp/owner.json').read_text())
        assert failed['status'] == 'terminal' and not failed['passed'] and not failed['remaining_owned_processes']
        paths += [W / 'machine-dp/owner.json', W / 'machine-dp/stderr.log',
                  W / 'target-feature-review.json', R / product['ir']['path']]
        assembler = json.loads((W / 'toolchain.json').read_text())['tools']['llvm-as']
        assert pin(R / assembler['resolved']) == {k: assembler[k] for k in ('bytes', 'sha256')}
        paths.append(R / assembler['resolved'])
    inputs = {str(p.relative_to(R)) if p.is_relative_to(R) else str(p): pin(p) for p in paths}
    argv = [str(D / 'llc'), '-O2', '-filetype=obj', '-relocation-model=pic',
            '-mcpu=skylake-avx512', *([] if features else ['-mattr=-fma']),
            str(B / 'dp-value.features.bc') if features else str(R / product['output']['path']),
            '-o', str(target)]
    write(plan_path, dict(argv=argv, inputs=inputs, repair_feature_attributes=features,
                         original_optimized_IR=product['ir']))
    actions = json.loads((W / 'actions.json').read_text())
    assert case not in actions and not target.exists()
    actions[case] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=['/usr/bin/python3', '-B', str(W / 'machine.py'), 'emit-features' if features else 'emit'],
        frozen_inputs=sorted(inputs) + [str(plan_path.relative_to(R))],
        scope='One native lowering of the admitted strict DP value-interface module, plus static whole-object disassembly and symbol inventory. No execution, speed, numerical/provider/full scientific or release admission.')
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=case, compiler_run=False)))
elif base_phase == 'emit':
    plan = json.loads(plan_path.read_text())
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    assert not target.exists()
    if features:
        original = plan['original_optimized_IR']
        assert pin(R / original['path']) == {k: original[k] for k in ('bytes', 'sha256')}
        text = (R / original['path']).read_text()
        token = ' "target-features"="-fma"'
        groups = text.count(token)
        assert groups > 0
        repaired = text.replace(token, '')
        # No instruction, FP metadata, mask, ABI, layout, count or source data
        # changes. Keep original CPU and all original C++ feature strings.
        assert len(text) - len(repaired) == groups * len(token)
        ll = B / 'dp-value.features.ll'
        bc = B / 'dp-value.features.bc'
        assert not ll.exists() and not bc.exists()
        with ll.open('x') as stream:
            stream.write(repaired)
        subprocess.run([str(D / 'llvm-as'), str(ll), '-o', str(bc)], check=True)
        write(W / 'target-repair.json', dict(original=original, removed_attribute_token=token,
            attribute_groups=groups, repaired_IR=dict(path=str(ll.relative_to(R)), **pin(ll)),
            bitcode=dict(path=str(bc.relative_to(R)), **pin(bc)), instruction_or_FP_changes=False,
            reason='Restore SKX AVX512 capability closure after -fma also disabled its dependent AVX512 features. Strict FP/no-contract semantics and actual no-FMA machine gate remain.'))
    subprocess.run(plan['argv'], check=True)
    listing = subprocess.check_output(['/usr/bin/objdump', '-d', '-M', 'intel', '--no-show-raw-insn', str(target)], text=True)
    ops = Counter(re.findall(r'^\s*[0-9a-f]+:\s+(\S+)', listing, re.M))
    assert not any(op.startswith(('vfm', 'vfnm')) or op in ('(bad)', '.byte') for op in ops), ops
    symbols = subprocess.check_output(['/usr/bin/nm', '-S', '--defined-only', str(target)], text=True)
    cores = [line for line in symbols.splitlines() if ' sirius_value_core_' in line]
    assert len(cores) == 8, cores
    undefined = subprocess.check_output(['/usr/bin/nm', '-u', str(target)], text=True).splitlines()
    result = dict(output=dict(path=str(target.relative_to(R)), **pin(target)), argv=plan['argv'],
        no_FMA=True, decoded_instructions=sum(ops.values()), instruction_counts=dict(ops),
        core_symbols=cores, undefined_symbols=undefined,
        representation_IR_gate=json.loads((W / 'value-optimize-joins-result.json').read_text())['output'],
        scope='Object generation and static inventory only. Aggregate ABI may use stack/spills. Exact primitive retention, full ISA closure, numerical equivalence, PC/FP-scope joins, provider and frame admission remain unestablished.')
    write(W / ('machine-features-result.json' if features else 'machine-result.json'), result)
    print(json.dumps({k: result[k] for k in ('output', 'no_FMA', 'decoded_instructions', 'core_symbols', 'undefined_symbols')}))
else:
    raise ValueError(phase)
