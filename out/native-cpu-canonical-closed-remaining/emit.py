"""Derive once, then emit Transport and Endpoint sequentially; no native run."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
V = R / 'out/native-cpu-strict-ir-linking'
GENERATOR = V / 'canonical_private_remaining.py'
CONTRACT = V / 'canonical_private_remaining_contract.json'
STAGES = ('transport', 'endpoint')
CONTRACT_PIN = dict(bytes=58311, sha256='f2f9b0ec0fb0930485eaf0016d219327044f65602d22a9696b47e9450f98e1cc')
FUTURE_UNWIND_GATE = ('Before the first native emission, the exact newly pinned strict/context result must carry '
    'asynchronous uwtable on the final ordinary private entry; verify its whole-function FDE/CFA/GP rules at all '
    'calls and returns before provider/C8. This ordinary-IR phase implements no strict/native/FDE pass.')


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def directory(stage):
    return R / ('bin/linux-gcc/native-cpu-canonical-closed-' + stage)


def source_path(stage):
    return directory(stage) / 'ispc-source' / (stage + '.ispc')


def terminal(name):
    path = W / name / 'owner.json'
    value = doc(path)
    assert value['status'] == 'terminal' and value['passed'] and value['child_exit'] == 0
    assert value['source_unchanged'] and value['whole_inputs_unchanged']
    assert not value['remaining_owned_processes'] and not value['cleanup_errors']
    assert value['observed_births'] and all(row['absent'] for row in value['observed_births'])
    assert path.with_name('inputs-before.json').read_bytes() == path.with_name('inputs-after.json').read_bytes()


def authority():
    assert pin(CONTRACT) == CONTRACT_PIN
    value = doc(CONTRACT)
    assert pin(GENERATOR) == value['generator']
    for name, expected in value['source_inputs'].items():
        assert pin(R / name) == expected, name
    for stage, tuples, layers, roots in (('transport', 3144, 159, 40), ('endpoint', 4148, 195, 100)):
        proof = value['source_proofs'][stage]
        for key in ('source_inverse_to_publication_exact', 'source_inverse_to_AS_exact',
            'unchanged_outside_export_binding_and_EvaluateProgram', 'generic_fallback_removed_only_from_new_private_body',
            'canonical_tuple_records_match_exact_tables', 'complete_context_memset_and_full8_physical_packing_unchanged',
            'supported_callers_lane0', 'mixed_families_supported', 'homogeneous_all_active_guard_absent',
            'all_root_validation_before_any_publication', 'original_combined_mask_five_word_store_order'):
            assert proof[key] is True, (stage, key)
        assert (proof['original_tuple_records'], proof['original_layer_call_sites'], proof['root_positions']) == (tuples, layers, roots)
        assert proof['publication_word_store_sites'] == roots * 5
        assert proof['per_layer_family_order'] == sorted(proof['per_layer_family_order'], key=lambda row: (row[0], row[1]))
        assert len(proof['per_layer_family_order']) == layers
        if stage == 'endpoint':
            for key in ('live_Endpoint11241_guard_unchanged', 'Endpoint_prefix_state_and_whole_state_rollback_unchanged',
                        'required20_and100_supported_per_lane'):
                assert proof[key] is True
            assert proof['required_roots'] == [20, 100]
        else:
            assert proof['required_roots'] == 40
        argv = value['prospective_emit_commands'][stage]
        assert argv[1] == str(source_path(stage))
        assert argv[-2:] == ['-o', str(directory(stage) / (stage + '.ll'))]
        assert argv[2:-2] == ['-O2', '--target=avx512skx-x8', '--cpu=skx', '--arch=x86-64',
            '--addressing=32', '--PIC', '--opt=disable-fma', '--math-lib=system', '--emit-llvm-text']
    return value


def generated_sources(contract):
    receipt = doc(W / 'derive-result.json')
    assert receipt['source_proofs'] == contract['source_proofs']
    for stage in STAGES:
        expected = dict(path=str(source_path(stage).relative_to(R)), **contract['derived_sources'][stage])
        assert receipt['outputs'][stage] == expected
        assert pin(source_path(stage)) == contract['derived_sources'][stage]


phase, = sys.argv[1:]
contract = authority()
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    assert all(not directory(stage).exists() for stage in STAGES)
    paths = [*W.glob('*.py'), GENERATOR, CONTRACT,
        R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
        R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py']
    inputs = dict(contract['source_inputs'])
    inputs.update({str(path.relative_to(R)): pin(path) for path in paths})
    commands = dict(derive=['/usr/bin/python3', '-B', str(GENERATOR), 'derive'],
        **{'emit-' + stage: contract['prospective_emit_commands'][stage] for stage in STAGES})
    write(W / 'plan.json', dict(commands=commands, inputs=inputs,
        source_contract=CONTRACT_PIN, stage_order=list(commands), future_native_unwind_gate=FUTURE_UNWIND_GATE))
    actions = {}
    for name in commands:
        extra = []
        if name != 'derive':
            extra += [str(source_path(stage).relative_to(R)) for stage in STAGES]
            extra += [str((W / path).relative_to(R)) for path in ('derive-result.json', 'derive/owner.json',
                'derive/inputs-before.json', 'derive/inputs-after.json')]
        if name == 'emit-endpoint':
            extra += [str((W / path).relative_to(R)) for path in ('transport-inventory.json',
                'emit-transport/owner.json', 'emit-transport/inputs-before.json', 'emit-transport/inputs-after.json')]
            extra += [str((directory('transport') / 'transport.ll').relative_to(R))]
        actions[name] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            kind='diagnostic', expected_exit=0, software_provider=False,
            limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
            argv=['/usr/bin/python3', '-B', str(Path(__file__).resolve()), name],
            frozen_inputs=sorted(inputs) + [str((W / 'plan.json').relative_to(R))] + extra,
            scope='Closed mixed canonical T/E source derivation and ordinary ISPC IR emission only, ordered derive→Transport→Endpoint and stop first failure. Original tables/tuples/mixed masks/refusals/Endpoint11241/FP/source/count/layout/governors remain. No strict/native/provider/runtime/ABBA/frame or scientific acceptance. Future native code must have checked async unwind metadata from its first emission.')
    write(W / 'actions.json', actions)
    print(json.dumps(dict(prepared=list(commands), runtime_run=False)))
else:
    plan = doc(W / 'plan.json')
    assert phase in plan['commands'] and plan['stage_order'] == ['derive', 'emit-transport', 'emit-endpoint']
    assert plan['source_contract'] == CONTRACT_PIN
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    if phase != 'derive':
        terminal('derive')
        generated_sources(contract)
        if phase == 'emit-endpoint':
            terminal('emit-transport')
            previous = doc(W / 'transport-inventory.json')['output']
            assert previous == dict(path=str((directory('transport') / 'transport.ll').relative_to(R)),
                **pin(directory('transport') / 'transport.ll'))
        stage = phase.removeprefix('emit-')
        assert stage in STAGES and not (directory(stage) / (stage + '.ll')).exists()
    subprocess.run(plan['commands'][phase], check=True)
    if phase == 'derive':
        outputs = {}
        for stage in STAGES:
            assert pin(source_path(stage)) == contract['derived_sources'][stage]
            outputs[stage] = dict(path=str(source_path(stage).relative_to(R)), **pin(source_path(stage)))
        write(W / 'derive-result.json', dict(outputs=outputs, source_proofs=contract['source_proofs'],
            source_contract=CONTRACT_PIN, compiler_run=False, numerical_qualification_transferred=False))
    else:
        path = directory(stage) / (stage + '.ll')
        text = path.read_text()
        definitions = re.findall(r'^define .*$', text, re.M)
        assert definitions
        ordinary = Counter(re.findall(r'^\s+%[^=]+ = (fadd|fsub|fmul|fdiv|frem|fneg|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b', text, re.M))
        opcodes = Counter(re.findall(r'^  (?:%[^=\n]+ = )?(?:(?:tail|musttail|notail)\s+)?([a-z][a-z0-9.]*)\b', text, re.M))
        calls = Counter(re.findall(r'\bcall\s+.*?@([^ (]+)\(', text))
        predicates = Counter(re.findall(r'= fcmp (\w+) ', text))
        fast = re.findall(r'^\s+.*\b(?:fadd|fsub|fmul|fdiv|frem|fneg|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b.*$', text, re.M)
        assert not fast, fast[:3]
        assert not any(name.startswith(('llvm.fma.', 'llvm.fmuladd.')) for name in calls)
        result = dict(output=dict(path=str(path.relative_to(R)), **pin(path)), stage=stage,
            source=dict(path=str(source_path(stage).relative_to(R)), **pin(source_path(stage))),
            definitions=len(definitions), definition_headers=definitions, instruction_opcodes=dict(opcodes),
            ordinary_fp_operations=dict(ordinary), fp_instruction_flags={}, comparison_predicates=dict(predicates),
            call_inventory=dict(calls), math_calls={name: count for name, count in calls.items() if name.startswith('sirius_native_v4_')},
            source_proof=contract['source_proofs'][stage], future_native_unwind_gate=FUTURE_UNWIND_GATE,
            scope='Actual first ordinary ISPC IR inventory only; no copied expected function/FP/math call totals. Exact current mixed source proofs remain bound. No strict/body/native/numerical/benefit/provider/frame/scientific/release qualification.')
        write(W / (stage + '-inventory.json'), result)
        print(json.dumps({key: result[key] for key in ('output', 'stage', 'definitions', 'ordinary_fp_operations', 'comparison_predicates', 'math_calls')}))
