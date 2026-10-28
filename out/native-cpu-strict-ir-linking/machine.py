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
if phase == 'prepare':
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
    inputs = {str(p.relative_to(R)) if p.is_relative_to(R) else str(p): pin(p) for p in paths}
    argv = [str(D / 'llc'), '-O2', '-filetype=obj', '-relocation-model=pic',
            '-mcpu=skylake-avx512', '-mattr=-fma', str(R / product['output']['path']),
            '-o', str(B / 'dp-value.o')]
    write(W / 'machine-plan.json', dict(argv=argv, inputs=inputs))
    actions = json.loads((W / 'actions.json').read_text())
    assert 'machine-dp' not in actions and not (B / 'dp-value.o').exists()
    actions['machine-dp'] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=['/usr/bin/python3', '-B', str(W / 'machine.py'), 'emit'],
        frozen_inputs=sorted(inputs) + ['out/native-cpu-strict-ir-linking/machine-plan.json'],
        scope='One native lowering of the admitted strict DP value-interface module, plus static whole-object disassembly and symbol inventory. No execution, speed, numerical/provider/full scientific or release admission.')
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared='machine-dp', compiler_run=False)))
elif phase == 'emit':
    plan = json.loads((W / 'machine-plan.json').read_text())
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    target = B / 'dp-value.o'
    assert not target.exists()
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
    write(W / 'machine-result.json', result)
    print(json.dumps({k: result[k] for k in ('output', 'no_FMA', 'decoded_instructions', 'core_symbols', 'undefined_symbols')}))
else:
    raise ValueError(phase)
