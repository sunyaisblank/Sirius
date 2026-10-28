"""Attribute-only asynchronous unwind repair; one assemble/native action."""
from pathlib import Path
import hashlib
import importlib.util
import json
import re
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-closed-dp-unwind'
OLD = R / 'out/native-cpu-canonical-closed-dp'
V = R / 'out/native-cpu-strict-ir-linking'
IR = R / 'bin/linux-gcc/native-cpu-canonical-closed-dp/dp.closed.ll'
OBJ = R / 'bin/linux-gcc/native-cpu-canonical-closed-dp/dp.closed.o'
IR_PIN = dict(bytes=43874147, sha256='6a86140bbb8fa7f179255df906522dae6ea1cc38d3cd9e44f02084a724c8fdd9')
OBJ_PIN = dict(bytes=2850896, sha256='ae66bec69645ec96f91391d9e2471b1d2759cb682458fbe104c98a90060e6bf2')


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def terminal(path):
    row = doc(path)
    assert row['status'] == 'terminal' and row['passed'] and row['child_exit'] == 0
    assert row['source_unchanged'] and row['whole_inputs_unchanged']
    assert not row['remaining_owned_processes'] and not row['cleanup_errors']
    assert row['observed_births'] and all(x['absent'] for x in row['observed_births'])
    assert path.with_name('inputs-before.json').read_bytes() == path.with_name('inputs-after.json').read_bytes()


def candidate():
    assert pin(IR) == IR_PIN and pin(OBJ) == OBJ_PIN
    original = IR.read_text()
    definitions = re.findall(r'^define .*$', original, re.M)
    assert len(definitions) == 1 and definitions[0].startswith('define hidden void @sirius_closed_dp_packet(')
    assert definitions[0].endswith('#0 {') and len(re.findall(r'#0\b', original)) == 2
    before, = re.findall(r'^attributes #0 = .*$', original, re.M)
    assert before.startswith('attributes #0 = { nounwind strictfp memory(readwrite) ')
    assert 'uwtable' not in before
    after = before.replace('nounwind strictfp ', 'nounwind strictfp uwtable ', 1)
    assert original.count(before) == 1 and after not in original
    transformed = original.replace(before, after, 1)
    assert transformed.replace(after, before, 1) == original
    return transformed, before, after


phase, = sys.argv[1:]
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    assert not B.exists()
    text, before, after = candidate()
    terminal(OLD / 'native/owner.json')
    native = doc(OLD / 'native-result.json')
    assert native['output'] == dict(path=str(OBJ.relative_to(R)), **OBJ_PIN)
    tools = doc(V / 'toolchain.json')['tools']
    assembler = R / tools['llvm-as']['path']
    assert pin(assembler.resolve()) == {k: tools['llvm-as'][k] for k in ('bytes', 'sha256')}
    commands = [
        [str(assembler), str(B / 'dp.unwind.ll'), '-o', str(B / 'dp.unwind.bc')],
        list(native['argv'])]
    commands[1][commands[1].index(str(R / 'bin/linux-gcc/native-cpu-canonical-closed-dp/dp.closed.bc'))] = str(B / 'dp.unwind.bc')
    commands[1][commands[1].index('-o') + 1] = str(B / 'dp.unwind.o')
    assert pin(Path(commands[1][0]).resolve()) == {k: tools['llc'][k] for k in ('bytes', 'sha256')}
    paths = [*W.glob('*.py'), IR, OBJ, OLD / 'native-result.json', OLD / 'native/owner.json',
        OLD / 'native/inputs-before.json', OLD / 'native/inputs-after.json', OLD / 'cleanup-result.json',
        V / 'toolchain.json', R / 'out/native-cpu-field-compatibility/elf_proof.py',
        R / 'out/native-cpu-accumulate-compatibility/elf_proof.py',
        R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py', R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py',
        assembler.resolve(), Path(commands[1][0]).resolve(), Path('/usr/bin/readelf').resolve(), Path('/usr/bin/objdump').resolve()]
    inputs = {str(p.relative_to(R)) if p.is_relative_to(R) else str(p): pin(p) for p in paths}
    inputs.update(doc(OLD / 'context-plan.json')['inputs'])
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    write(W / 'build-plan.json', dict(inputs=inputs, commands=commands, original_native_argv=native['argv'],
        original_IR=dict(path=str(IR.relative_to(R)), **IR_PIN), original_object=dict(path=str(OBJ.relative_to(R)), **OBJ_PIN),
        attribute_before=before, attribute_after=after,
        derived_IR=dict(path=str((B / 'dp.unwind.ll').relative_to(R)), bytes=len(text.encode()), sha256=hashlib.sha256(text.encode()).hexdigest()),
        full_text_inverse=True, sole_attribute_group_user=True, UWTableKind='Async=Default=2',
        parser_source='https://raw.githubusercontent.com/llvm/llvm-project/85ac560262434c9ccfc0c183ec22d4138ed647fb/llvm/lib/AsmParser/LLParser.cpp#L2464-L2477',
        scope='Only sole DP entry uwtable attribute changes. Same name/CABI/IR instructions/FP metadata/masks/CFG and original llc flags. No original input/product is overwritten; no numerical acceptance transfer.'))
    write(W / 'actions.json', {'build': dict(
        revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=['/usr/bin/python3', '-B', str(W / 'build.py'), 'build'],
        frozen_inputs=sorted(inputs) + ['out/native-cpu-closed-dp-unwind/build-plan.json'],
        scope='One llvm-as23 and unchanged llc invocation for attribute-only async unwind tables, then read-only whole-function/relocation/FDE audit. No new math/raw optimization, linking, provider/workload or numerical/scientific/release qualification.')})
    print(json.dumps(dict(prepared=True, commands=2, only_function_attribute_changed=True)))
elif phase == 'build':
    plan = doc(W / 'build-plan.json')
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    text, before, after = candidate()
    assert (before, after) == (plan['attribute_before'], plan['attribute_after'])
    assert dict(bytes=len(text.encode()), sha256=hashlib.sha256(text.encode()).hexdigest()) == {k: plan['derived_IR'][k] for k in ('bytes', 'sha256')}
    B.mkdir()
    with (B / 'dp.unwind.ll').open('x') as stream:
        stream.write(text)
    for command in plan['commands']:
        subprocess.run(command, check=True)
    spec = importlib.util.spec_from_file_location('owned_unwind_verifier', W / 'verify.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    proof = module.verify(OBJ, B / 'dp.unwind.o')
    proof['unwind_table'].update(asynchronous=True, asynchronous_evidence='Bare uwtable maps to UWTableKind::Default=Async in pinned LLVM23.1.2 parser/header; exact source inverse bound by build plan.')
    write(W / 'code-proof.json', proof)
    write(W / 'build-result.json', dict(passed=True, output=dict(path=str((B / 'dp.unwind.o').relative_to(R)), **pin(B / 'dp.unwind.o')),
        bitcode=dict(path=str((B / 'dp.unwind.bc').relative_to(R)), **pin(B / 'dp.unwind.bc')),
        source_original_object=plan['original_object'], source_original_IR=plan['original_IR'],
        only_function_attribute_changed=True, full_text_inverse=True, compile_commands=plan['commands'],
        code_proof=dict(path=str((W / 'code-proof.json').relative_to(R)), **pin(W / 'code-proof.json')),
        numerical_qualification_transferred=False,
        scope='Attribute-only async unwind object with full original function/relocation identity and checked FDE/CFI. No provider/runtime/numerical qualification; failed original C8 remains separate evidence.'))
    print(json.dumps(dict(passed=True, only_function_attribute_changed=True, full_function_bytes_identical=True,
        normalized_code_relocations_identical=True, FDE_verified=True)))
else:
    raise ValueError(phase)
