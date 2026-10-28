"""Bind the compiled shared interface; transfer no numerical qualification."""
from pathlib import Path
import hashlib
import json
import re

R = Path.cwd()
C = R / 'out/native-cpu-shared-value-component'
V = R / 'out/native-cpu-strict-ir-linking'
OLD = R / 'out/native-cpu-dispatch-token-provider/provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
HEADERS = ('src/gallium/drivers/llvmpipe/lp_native_cpu_bridge.h',
           'src/gallium/drivers/llvmpipe/sirius_native_bridge.h')
CORES = {f'sirius_shared_value_{kind}_{op}8'
         for kind, ops in (('dense', ('add', 'multiply', 'divide', 'sqrt')),
                          ('tagged', ('add', 'subtract', 'multiply', 'divide', 'sqrt')))
         for op in ops}
INTERFACE = dict(compiled_interface_passed=True, static_symbol_binding_passed=True,
                 numerical_qualification_transferred=False, genuine_job_numerical_gate_required=True)


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def checked(path, expected):
    assert pin(path) == {k: expected[k] for k in ('bytes', 'sha256')}, path


def terminal(path):
    value = doc(path)
    assert value['status'] == 'terminal' and value['passed'] and value['child_exit'] == 0
    assert value['source_unchanged'] and value['whole_inputs_unchanged']
    assert not value['remaining_owned_processes'] and not value['cleanup_errors']
    assert value['observed_births'] and all(row['absent'] for row in value['observed_births'])
    assert path.with_name('inputs-before.json').read_bytes() == path.with_name('inputs-after.json').read_bytes()
    return value


def contract_inputs(value):
    for name, expected in value['source_inputs'].items():
        checked(R / name, expected)


def component():
    old = doc(OLD / 'source-contract.json')
    assert pin(OLD / 'source-contract.json') == dict(bytes=134818,
        sha256='16d09bcdbc1d8e853b8f5ac9666af890f24ee501989314622121a202214f7125')
    contract_inputs(old)
    headers = {name: old['modified_files'][name]['candidate'] for name in HEADERS}
    for name, expected in headers.items():
        checked(S / name, expected)
    terminal(C / 'build/owner.json')
    plan, result, derivation = (doc(C / name) for name in ('build-plan.json', 'build-result.json', 'derivation.json'))
    assert result['passed'] and result['all_raw_scope_and_core_references_bound']
    assert set(result['shared_core_definitions']) == CORES and len(result['shared_core_definitions']) == 9
    assert not any(name in result['remaining_external_symbols'] for name in ('sirius_shared_value_', 'sirius_scope_'))
    for name, expected in plan['inputs'].items():
        checked(R / name, expected)
    for name, expected in derivation['source_inputs'].items():
        checked(R / name, expected)
    checked(R / derivation['derived']['path'], derivation['derived'])
    assert derivation['public_wrapper_prefix_byte_preserved']
    # Recreate exactly the three callback substitutions and their declarations.
    original = (R / 'out/native-cpu-dispatch-admission-token/rows_adapter.cpp').read_text()
    boundary = original.index('#include "admitted_bridge.h"')
    public, private, declarations = original[:boundary], original[boundary:], []
    for stage, scope in (('transport', 'transport'), ('endpoint', 'endpoint'), ('dopri_phase', 'dp')):
        before, after = f'sirius_simd_v4_{stage}_packet(', f'sirius_scope_{scope}_packet('
        assert private.count(before) == 1
        private = private.replace(before, after, 1)
        declarations.append('extern "C" __attribute__((visibility("hidden"))) void ' + after +
            'uint32_t *input, uint32_t input_words, uint32_t *output, uint32_t output_words, uint32_t first_x, uint32_t groups);')
    for name in ('admitted_bridge.h', 'native_packet.hpp'):
        before = f'#include "{name}"'
        assert private.count(before) == 1
        private = private.replace(before, f'#include "../native-cpu-dispatch-admission-token/{name}"', 1)
    assert (C / 'rows_adapter.cpp').read_text() == public + '\n'.join(declarations) + '\n' + private
    adapter = dict(path='bin/linux-gcc/native-cpu-shared-value-component/adapter-only.o', **result['adapter'])
    checked(R / adapter['path'], adapter)
    expected_compile = list(derivation['compile_argv'])
    expected_compile[expected_compile.index('-o') + 1] = str(R / adapter['path'])
    assert plan['compile_argv'] == expected_compile
    reused = doc(R / 'out/native-cpu-dispatch-token-provider/merge-inputs.json')['parts'][1:]
    assert result['reused'] == plan['reused'] == reused and len(reused) == 6
    products, receipt_paths = [], []
    for stage in ('dp', 'transport', 'endpoint', 'math'):
        case = 'shared-endpoint-cleanup-native' if stage == 'endpoint' else f'shared-{stage}-native'
        owner, receipt = V / case / 'owner.json', V / (case + '-result.json')
        terminal(owner)
        value = doc(receipt)
        assert value['no_FMA'] and value['decoded_instructions'] > 0
        checked(R / value['output']['path'], value['output'])
        products.append(value['output'])
        receipt_paths += [owner, owner.with_name('inputs-before.json'), owner.with_name('inputs-after.json'), receipt]
    assert plan['new_value_objects'] == products
    assert [row['object'] for row in result['raw_object_references']] == products
    refs = [set(row['shared_core_references']) for row in result['raw_object_references']]
    assert refs == [CORES - {'sirius_shared_value_tagged_sqrt8'}, CORES, CORES, set()]
    abi_review = doc(C / 'abi-review.json')
    assert set(abi_review['objects']) == {'shared-dp-native', 'shared-transport-native',
        'shared-endpoint-cleanup-native', 'shared-math-native'}
    assert set(row['sha256'] for row in abi_review['objects'].values()) == set(row['sha256'] for row in products)
    assert abi_review['total_shared_core_calls'] == 8241
    assert abi_review['per_raw'] == dict(dp=1581, transport=2701, endpoint=3959)
    assert abi_review['no_false_formal_ABI_constraints']
    assert abi_review['public_wrapper_comparison']['nonrelocation_bytes_identical']
    assert abi_review['public_wrapper_comparison']['offsets_types_addends_target_names_identical']
    shared_plan = doc(V / 'shared-plan.json')
    assert set(shared_plan['core_names']) == CORES
    for name, expected in shared_plan['inputs'].items():
        checked(R / name, expected)
    math_ir = doc(V / 'shared-math-opt-result.json')
    signatures = math_ir['fixed_value_signatures']
    assert set(signatures) == CORES and len(math_ir['core_definitions']) == 9
    for header in math_ir['core_definitions']:
        attributes = math_ir['definition_attributes'][re.search(r'#(\d+)', header)[1]]
        assert 'noinline' in attributes and 'strictfp' in attributes
        effects = re.search(r'\bmemory\(([^)]+)\)', attributes)[1]
        assert effects in ('inaccessiblemem: readwrite',
            'read, argmem: none, inaccessiblemem: readwrite, errnomem: none, target_mem: none')
    core_calls = dict(dopri_phase=1581, transport=2701, endpoint=3959)
    structural_paths = [V / 'shared-plan.json', V / 'shared-math-opt-result.json']
    for stage, count in (('dp', 1581), ('transport', 2701), ('endpoint-cleanup', 3959), ('math', None)):
        case = f'shared-{stage}-opt'
        receipt = V / (case + '-result.json')
        terminal(V / case / 'owner.json')
        ir = doc(receipt)
        assert ir['payload_allocas'] == ir['tag_allocas'] == ir['adapter_calls'] == 0
        assert ir['constrained_calls']
        if count is not None:
            assert not ir['core_definitions'] and sum(ir['core_calls'].values()) == count
        for header in ir['core_declarations']:
            attributes = ir['definition_attributes'][re.search(r'#(\d+)', header)[1]]
            assert 'noinline' in attributes and 'strictfp' in attributes
            assert re.search(r'\bmemory\(([^)]+)\)', attributes)[1] == (
                'read, argmem: none, inaccessiblemem: readwrite, errnomem: none, target_mem: none')
        for name, signature in ir['fixed_value_signatures'].items():
            assert name in signatures
            assert {k: v for k, v in signature.items() if k != 'header'} == {
                k: v for k, v in signatures[name].items() if k != 'header'}
        for row in (ir['output'], ir['ir']):
            checked(R / row['path'], row)
            structural_paths.append(R / row['path'])
        structural_paths += [receipt, V / case / 'owner.json', V / case / 'inputs-before.json', V / case / 'inputs-after.json']
    merged = result['output']
    assert merged['path'] == 'bin/linux-gcc/native-cpu-shared-value-component/rows_adapter.o'
    parts = [adapter] + reused + products
    assert len(parts) == 11
    assert plan['merge_argv'] == ['/usr/bin/ld', '-r', '-o', str(R / merged['path'])] + [str(R / row['path']) for row in parts]
    for row in parts + [merged] + old['native_objects'][:7]:
        checked(R / row['path'], row)
    paths = [C / name for name in ('build-plan.json', 'build-result.json', 'derivation.json', 'abi-review.json', 'build/owner.json',
                                  'build/inputs-before.json', 'build/inputs-after.json')]
    paths += [OLD / 'source-contract.json', R / 'out/native-cpu-dispatch-token-provider/merge-inputs.json'] + receipt_paths
    inputs = dict(plan['inputs'])
    inputs.update(shared_plan['inputs'])
    paths += structural_paths
    inputs.update({str(path.relative_to(R)): pin(path) for path in paths})
    inputs.update({str((S / name).relative_to(R)): expected for name, expected in headers.items()})
    return dict(native_first7=old['native_objects'][:7], merged=merged, parts=parts,
        raw_v4_objects=[products[1], products[2], products[0]], shared_value_math_object=products[3],
        shared_value_interface=dict(passed=True, core_count=9, raw_core_calls=core_calls,
            fixed_signatures_preserved=True, no_ordinary_writes=True, strict_scope_source_bound=True),
        component_interface=dict(INTERFACE), unchanged_bridge_headers=headers, source_inputs=inputs,
        scope='Compiled adapter/shared-value interface and static symbol binding only. No complete linked-PC/ISA/FP-scope or numerical acceptance transfers. First numerical gate requires current genuine-positive-job C8 with the full original FP scope; then sole ABBA/Decimal, all original controls and unchanged frame. Direct private fixture calls outside that job scope are not new numerical qualification.')
