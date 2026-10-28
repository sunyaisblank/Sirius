"""Bind mixed AU/closed-DP compiled interfaces; transfer no numerical pass."""
from pathlib import Path
import importlib.util
import re

R = Path.cwd()
C = R / 'out/native-cpu-closed-dp-component'
AU = R / 'out/native-cpu-shared-value-provider'
DP = R / 'out/native-cpu-canonical-closed-dp'
M = R / 'out/native-cpu-local-effects'
V = R / 'out/native-cpu-strict-ir-linking'
_path = AU / 'component_gate.py'
import hashlib
assert hashlib.sha256(_path.read_bytes()).hexdigest() == '197dc4b7e74298af216414daa72ee7efb0d14d42b896a0299cbee31050789d7c'
_spec = importlib.util.spec_from_file_location('AV_pinned_AU_component', _path)
_au = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_au)
doc, pin, checked, terminal, contract_inputs = _au.doc, _au.pin, _au.checked, _au.terminal, _au.contract_inputs
HEADERS, CORES = _au.HEADERS, _au.CORES
LOCAL_CORES = {name.replace('sirius_shared_value_', 'sirius_local_value_') for name in CORES}
DP_CALLS = dict(sirius_local_value_dense_add8=1, sirius_local_value_tagged_add8=585,
    sirius_local_value_tagged_subtract8=284, sirius_local_value_tagged_multiply8=661,
    sirius_local_value_tagged_divide8=45)


def local_abi(math_text, raw_text):
    """Exact integer-only cc92 signatures and strict inaccessible effects."""
    def headers(text, defined):
        attrs = dict(re.findall(r'^attributes #(\d+) = (.*)$', text, re.M))
        rows = {}
        for header in re.findall(r'^' + ('define' if defined else 'declare') + r' .*@sirius_local_value_.*$', text, re.M):
            match = re.search(r'@([^ (]+)\(([^)]*)\)', header)
            name, args = match.groups()
            assert name in LOCAL_CORES and name not in rows
            tagged = '_tagged_' in name
            count = 12 if tagged else 10
            ret = '{ ' + ', '.join(['<8 x i32>'] * (6 if tagged else 5)) + ' }'
            assert header.startswith(('define' if defined else 'declare') + ' hidden x86_regcallcc ' + ret + ' @' + name + '(')
            clean_args = re.sub(r' %[-a-zA-Z$._0-9]+', '', args)
            assert clean_args == ', '.join(['<8 x i32>'] * count + ['i32']), header
            effect = attrs[re.search(r'#(\d+)', header)[1]]
            assert 'strictfp' in effect and 'noinline' in effect
            assert re.search(r'\bmemory\(([^)]+)\)', effect)[1] == 'inaccessiblemem: readwrite'
            rows[name] = dict(arguments=count, returns=6 if tagged else 5, calling_convention=92,
                integer_vectors_only=True, no_pointer_or_false_formal_constraints=True)
        return rows
    definitions, declarations = headers(math_text, True), headers(raw_text, False)
    assert set(definitions) == LOCAL_CORES and set(declarations) == set(DP_CALLS)
    for name, row in declarations.items():
        assert row == definitions[name]
    for pattern in (r'^target datalayout = .*$', r'^target triple = .*$'):
        assert re.search(pattern, math_text, re.M)[0] == re.search(pattern, raw_text, re.M)[0]
    assert 'metadata !"round.dynamic"' in math_text and 'metadata !"fpexcept.strict"' in math_text
    assert not re.search(r'metadata !"(?:round.tonearest|fpexcept.ignore|fpexcept.maytrap)"', math_text + raw_text)
    return definitions


def component():
    # AU authorizes only byte-reused constituents, not this new DP path.
    au = _au.component()
    terminal(C / 'build-audit/owner.json')
    plan, result, derivation = (doc(C / name) for name in ('build-plan.json', 'build-result.json', 'derivation.json'))
    audit = doc(C / 'build-audit-plan.json')
    failed_owner = dict(bytes=1840, sha256='33450b323f08bdea73d0ef3348f28798c19e690f674488728713adba26d1463b')
    original_plan = dict(bytes=21530, sha256='dc19ffb3ab77ea20cea39a8dfabb806e7b69355b5a025f07db823a62714cbaa8')
    retained = {
        'bin/linux-gcc/native-cpu-closed-dp-component/adapter-only.o': dict(bytes=7680,
            sha256='884c89b1bc837009efa855cf4eb198c0f17caf622e3021ee4097cfeb9ed3bdf6'),
        'bin/linux-gcc/native-cpu-closed-dp-component/rows_adapter.o': dict(bytes=22724920,
            sha256='9ffdb4519cb8065b49b0b25b8ae539bef24d300fe14becfd86b55368c27bdbb7')}
    checked(C / 'build/owner.json', failed_owner)
    checked(C / 'build-plan.json', original_plan)
    failed = doc(C / 'build/owner.json')
    assert failed['status'] == 'terminal' and failed['passed'] is False and failed['child_exit'] == 1
    assert failed['source_unchanged'] and failed['whole_inputs_unchanged']
    assert not failed['remaining_owned_processes'] and not failed['cleanup_errors']
    assert failed['observed_births'] and all(row['absent'] for row in failed['observed_births'])
    assert (C / 'build/inputs-before.json').read_bytes() == (C / 'build/inputs-after.json').read_bytes()
    checked(C / 'build/stderr.log', dict(bytes=400,
        sha256='ebbd36822d58b1a055b2c5ddd1189f30222f0265d9d06267aa66f07163f03c08'))
    assert audit['original_failed_owner'] == result['original_failed_owner'] == failed_owner
    assert audit['original_build_plan'] == result['original_build_plan'] == original_plan
    assert audit['original_compile_argv'] == plan['compile_argv'] and audit['original_merge_argv'] == plan['merge_argv']
    assert audit['retained_products'] == retained
    assert audit['accepted_artifact_owner'] == result['accepted_artifact_owner'] == 'out/native-cpu-closed-dp-component/build-audit/owner.json'
    assert not audit['compilation_repeated'] and not result['compilation_repeated']
    assert not audit['relocation_link_repeated'] and not result['relocation_link_repeated']
    checked(C / 'build-audit-plan.json', result['audit_plan'])
    reader = 'out/native-cpu-closed-dp-component/build.py'
    assert audit['source_hash_override'] == dict(path=reader, before=plan['inputs'][reader], after=pin(R / reader))
    for name, expected in plan['inputs'].items():
        assert audit['inputs'][name] == (pin(R / reader) if name == reader else expected), name
    for name, expected in audit['inputs'].items():
        checked(R / name, expected)
    for name, expected in retained.items():
        checked(R / name, expected)
    original_actions = doc(C / 'build/inputs-before.json')['out/native-cpu-closed-dp-component/actions.json']
    checked(C / 'build-audit-prior-actions.json', original_actions)
    assert doc(C / 'actions.json')['build'] == doc(C / 'build-audit-prior-actions.json')['build']
    assert result['passed'] and result['compiled_interface_passed'] and not result['numerical_qualification_transferred']
    assert result['all_raw_scope_and_core_references_bound']
    assert set(result['core_definitions']) == CORES | LOCAL_CORES and len(result['core_definitions']) == 18
    assert result['closed_raw_exports'] == ['sirius_closed_dp_packet']
    contract_inputs(derivation)
    checked(R / derivation['derived']['path'], derivation['derived'])
    assert derivation['source_inverse_exact'] and derivation['public_wrapper_prefix_byte_preserved']
    old_source = (R / 'out/native-cpu-shared-value-component/rows_adapter.cpp').read_text()
    candidate = (C / 'rows_adapter.cpp').read_text()
    assert old_source.count('sirius_scope_dp_packet(') == candidate.count('sirius_closed_dp_packet(') == 2
    assert candidate.replace('sirius_closed_dp_packet(', 'sirius_scope_dp_packet(') == old_source
    adapter = dict(path='bin/linux-gcc/native-cpu-closed-dp-component/adapter-only.o', **result['adapter'])
    checked(R / adapter['path'], adapter)
    assert result['adapter'] == retained[adapter['path']]
    assert plan['compile_argv'] == derivation['compile_argv']
    reused, au_values = au['parts'][1:7], [*au['raw_v4_objects'][:2], au['shared_value_math_object']]
    assert result['reused'] == plan['reused'] == reused and len(reused) == 6
    assert result['reused_au_value_objects'] == plan['reused_au_value_objects'] == au_values
    paths = []
    for root, stages in ((DP, ('context', 'cleanup', 'native')), (M, ('math-lower', 'math-native'))):
        for stage in stages:
            terminal(root / stage / 'owner.json')
            receipt = root / (stage + '-result.json')
            value = doc(receipt)
            checked(R / value['output']['path'], value['output'])
            paths += [receipt, root / stage / 'owner.json', root / stage / 'inputs-before.json', root / stage / 'inputs-after.json']
    native, cleanup = doc(DP / 'native-result.json'), doc(DP / 'cleanup-result.json')
    math_native, math_ir = doc(M / 'math-native-result.json'), doc(M / 'math-lower-result.json')
    assert native['output'] == dict(path='bin/linux-gcc/native-cpu-canonical-closed-dp/dp.closed.o',
        bytes=2850896, sha256='ae66bec69645ec96f91391d9e2471b1d2759cb682458fbe104c98a90060e6bf2')
    assert math_native['output'] == dict(path='bin/linux-gcc/native-cpu-local-effects/math.o',
        bytes=259216, sha256='69f9d771c5cf2d0de63e32ac70ac534961e9283a02753047784902e3ee54cdbe')
    assert native['core_relocations'] == cleanup['value_calls'] == DP_CALLS
    assert native['instructions'] == 321681 and not native['emitted_fused_instructions']
    assert cleanup['definitions'] == 1 and cleanup['owned_context_allocas'] == 0
    assert cleanup['memory_instruction_counts'] == dict(alloca=1, load=5)
    assert cleanup['masked_store_calls'] == 1010
    assert not math_native['emitted_fused_instructions'] and math_ir['ordinary_FP_forbidden']
    for row in (cleanup['ir'], math_ir['ir']):
        checked(R / row['path'], row)
        paths.append(R / row['path'])
    raw_text, math_text = (R / cleanup['ir']['path']).read_text(), (R / math_ir['ir']['path']).read_text()
    definitions = re.findall(r'^define .*$', raw_text, re.M)
    assert len(definitions) == 1 and definitions[0].startswith('define hidden void @sirius_closed_dp_packet(')
    assert '___un_' not in definitions[0]
    params = re.search(r'@sirius_closed_dp_packet\((.*)\) local_unnamed_addr', definitions[0])[1].split(', ')
    assert len(params) == 6
    assert all(row.startswith('ptr ' if index in (0, 2) else 'i32 ') for index, row in enumerate(params))
    signatures = local_abi(math_text, raw_text)
    memory_proof_path = V / 'private_finite_closure.json'
    memory_proof = doc(memory_proof_path)
    assert memory_proof['decision'] == 'GO_REMOVE_ONLY_OTHER_READ'
    assert not memory_proof['finite_FP_obligation']['ordinary_FP_conversion_authorized']
    assert not memory_proof['finite_FP_obligation']['memory_none_authorized']
    for name, expected in memory_proof['source_inputs'].items():
        checked(R / name, expected)
    context_contract = doc(V / 'closed_dp_context_contract.json')
    for name, expected in context_contract['source_inputs'].items():
        checked(R / name, expected)
    assert context_contract['phase'] == 'context' and context_contract['input']['sha256'] == '3f69125b95e0a6282a3e7dfa2d7617a0ea933af22dd2aa411314a8b1c76a2cf0'
    checked(R / context_contract['output']['path'], context_contract['output'])
    source_proof = doc(V / 'canonical_private_dp_contract.json')['source_proof']
    assert source_proof['sole_supported_caller_required200_lane0'] and source_proof['all_validation_before_first_publication']
    assert source_proof['complete_private_context_memset_preserved'] and source_proof['full8_physical_payloads_and_tags_unchanged']
    added = [native['output'], math_native['output']]
    assert result['new_closed_objects'] == plan['new_closed_objects'] == added
    merged = result['output']
    assert merged['path'] == 'bin/linux-gcc/native-cpu-closed-dp-component/rows_adapter.o'
    assert {key: merged[key] for key in ('bytes', 'sha256')} == retained[merged['path']]
    parts = [adapter] + reused + au_values + added
    assert len(parts) == 12
    assert plan['merge_argv'] == ['/usr/bin/ld', '-r', '-o', str(R / merged['path'])] + [str(R / row['path']) for row in parts]
    for row in parts + [merged]:
        checked(R / row['path'], row)
    inputs = dict(au['source_inputs']) | dict(audit['inputs']) | dict(context_contract['source_inputs'])
    paths += [C / name for name in ('build-plan.json', 'build-result.json', 'derivation.json', 'build/owner.json', 'build/inputs-before.json', 'build/inputs-after.json',
        'build-audit-plan.json', 'build-audit-prior-actions.json', 'build-audit/owner.json', 'build-audit/inputs-before.json', 'build-audit/inputs-after.json')]
    paths += [Path(__file__).resolve(), _path, V / 'closed_dp_context_contract.json', V / 'canonical_private_dp_contract.json',
        R / context_contract['output']['path'], memory_proof_path, AU / 'provider/source-contract.json', AU / 'provider/linked-native-proof.json']
    inputs.update({str(path.relative_to(R)): pin(path) for path in paths})
    interface = dict(compiled_interface_passed=True, static_symbol_binding_passed=True,
        numerical_qualification_transferred=False, genuine_job_numerical_gate_required=True)
    closed_interface = dict(passed=True, core_count=9, raw_core_calls=dict(dopri_phase=1576), fixed_signatures=signatures,
        fixed_signatures_preserved=True, no_ordinary_writes=True, no_FP_observers=True,
        strict_scope_source_bound=True, context_allocas=0, ordinary_loads=5, ordinary_stores=0,
        private_masked_stores=1010, sole_raw_export='sirius_closed_dp_packet', complete_source_initialization_proof_bound=True,
        actual_positive_job_required=True, complete_outer_FP_restore_required=True)
    return dict(native_first7=au['native_first7'], merged=merged, parts=parts,
        raw_v4_objects=au['raw_v4_objects'][:2] + [added[0]], shared_value_math_object=au['shared_value_math_object'],
        closed_dp_math_object=added[1], closed_dp_interface=closed_interface,
        shared_value_interface=dict(au['shared_value_interface'], raw_core_calls=dict(transport=2701, endpoint=3959)),
        component_interface=interface, unchanged_bridge_headers=au['unchanged_bridge_headers'], source_inputs=inputs,
        scope='AV compiled-interface/static source binding only. AU scoped T/E+shared9 and all legacy/public constituents retained exactly; admitted DP v4 alone selects closed DP with existing local-effects cc92 math9. No numerical transfer. Next gate is original genuine-job C8/all32; no ABBA/frame until remaining-stage closure is justified.')
