"""Bind the exact AV constituents plus one attribute-only DP unwind repair."""
from pathlib import Path
import hashlib
import importlib.util

R = Path.cwd()
C = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-closed-dp-unwind-component'
UW = R / 'out/native-cpu-closed-dp-unwind'
AV = R / 'out/native-cpu-closed-dp-provider'
_path = AV / 'component_gate.py'
assert hashlib.sha256(_path.read_bytes()).hexdigest() == 'db936d7b10703758849a8e43153f5e53747fb980211857b609c5bb9f2a9f4155'
_spec = importlib.util.spec_from_file_location('AW_pinned_AV_component', _path)
_av = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_av)
doc, pin, checked, terminal, contract_inputs = _av.doc, _av.pin, _av.checked, _av.terminal, _av.contract_inputs
HEADERS, CORES, LOCAL_CORES, DP_CALLS = _av.HEADERS, _av.CORES, _av.LOCAL_CORES, _av.DP_CALLS


def repaired():
    old = _av.component()
    terminal(UW / 'build/owner.json')
    result = doc(UW / 'build-result.json')
    assert result['passed'] is True
    proof_path = UW / 'code-proof.json'
    checked(proof_path, result['code_proof'])
    proof = doc(proof_path)
    assert result['only_function_attribute_changed'] is True
    assert result['numerical_qualification_transferred'] is False
    for key in ('passed', 'full_function_bytes_identical', 'normalized_code_relocations_identical',
                'sole_hidden_six_argument_C_ABI_preserved'):
        assert proof[key] is True, key
    table = proof['unwind_table']
    for key in ('full_function_FDE', 'asynchronous', 'complete_GP_save_rules', 'shrinkwrap_and_epilogues_checked'):
        assert table[key] is True, key
    previous = old['raw_v4_objects'][2]
    assert result['source_original_object'] == previous
    checked(R / previous['path'], previous)
    checked(R / result['source_original_IR']['path'], result['source_original_IR'])
    raw = result['output']
    assert raw['path'] == 'bin/linux-gcc/native-cpu-closed-dp-unwind/dp.unwind.o'
    checked(R / raw['path'], raw)
    parts = list(old['parts'])
    assert len(parts) == 12 and parts.index(previous) == 10
    assert parts[0]['path'] == 'bin/linux-gcc/native-cpu-closed-dp-component/adapter-only.o'
    assert parts[0]['bytes'] == 7680 and parts[0]['sha256'] == '884c89b1bc837009efa855cf4eb198c0f17caf622e3021ee4097cfeb9ed3bdf6'
    parts[10] = raw
    for row in parts:
        checked(R / row['path'], row)
    plan = doc(UW / 'build-plan.json')
    inputs = dict(old['source_inputs'])
    # Unwind production inputs are immutable compiler/source authorities;
    # Mesa's three identity/path targets are not inputs to this object repair.
    for name, expected in plan['inputs'].items():
        assert not name.endswith(('src/gallium/targets/lavapipe/meson.build',
            'src/gallium/drivers/llvmpipe/lp_screen.c', 'src/gallium/frontends/lavapipe/lvp_device.c'))
        checked(R / name, expected)
        inputs[name] = expected
    paths = [UW / name for name in ('build-plan.json', 'build-result.json', 'code-proof.json',
        'build/owner.json', 'build/inputs-before.json', 'build/inputs-after.json')]
    paths += [path for path in UW.iterdir() if path.is_file() and path.suffix == '.py']
    paths += [Path(__file__).resolve(), _path]
    for path in paths:
        inputs[str(path.relative_to(R))] = pin(path)
    for row in parts + [previous, result['source_original_IR']]:
        inputs[row['path']] = {key: row[key] for key in ('bytes', 'sha256')}
    repair = dict(compiled_attribute_repair_passed=True, executable_bytes_and_code_relocations_identical=True,
        object=raw, previous_object=previous, code_proof=dict(path=str(proof_path.relative_to(R)), **pin(proof_path)),
        object_unwind_table=table, only_function_attribute_changed=True,
        numerical_qualification_transferred=False, current_original_C8_required=True)
    return dict(parent=old, parts=parts, raw=raw, inputs=inputs, unwind_repair=repair)


def component():
    value = repaired()
    terminal(C / 'build/owner.json')
    plan, result = doc(C / 'build-plan.json'), doc(C / 'build-result.json')
    assert plan['parts'] == value['parts'] == result['parts'] and len(plan['parts']) == 12
    assert plan['new_native_compilations'] == result['new_native_compilations'] == 0
    assert result['passed'] and result['compiled_interface_passed'] and result['all_raw_scope_and_core_references_bound']
    assert result['numerical_qualification_transferred'] is False
    assert result['adapter'] == value['parts'][0]
    assert set(result['core_definitions']) == CORES | LOCAL_CORES and len(result['core_definitions']) == 18
    assert result['closed_raw_exports'] == ['sirius_closed_dp_packet']
    assert plan['merge_argv'] == ['/usr/bin/ld', '-r', '-o', str(B / 'rows_adapter.o')] + [str(R / row['path']) for row in value['parts']]
    merged = result['output']
    assert merged['path'] == str((B / 'rows_adapter.o').relative_to(R))
    checked(R / merged['path'], merged)
    for name, expected in plan['inputs'].items():
        checked(R / name, expected)
    inputs = dict(value['inputs'])
    inputs.update({str((C / name).relative_to(R)): pin(C / name) for name in (
        'build-plan.json', 'build-result.json', 'build/owner.json', 'build/inputs-before.json', 'build/inputs-after.json')})
    old = value['parent']
    answer = dict(old, merged=merged, parts=value['parts'],
        raw_v4_objects=old['raw_v4_objects'][:2] + [value['raw']], source_inputs=inputs,
        dp_unwind_repair=value['unwind_repair'],
        scope='AW compiled-interface/unwind-metadata repair only. Every AV constituent and adapter is reused byte-exact except DP.o, whose complete function bytes/code relocations are independently identical and whose new asynchronous full-function FDE is checked. No C++/raw/math recompilation or numerical transfer; original current genuine-job C8/all32 remains required.')
    return answer
