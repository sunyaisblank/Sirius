"""AV static verifier: pinned AU joins with an explicit mixed DP interface."""
from pathlib import Path
import hashlib

R = Path.cwd()
W = Path(__file__).resolve().parent
authority = R / 'out/native-cpu-shared-value-provider/verify_provider.py'
assert hashlib.sha256(authority.read_bytes()).hexdigest() == '55bf38021c9b70292338e9470e4105d7910114c1f8b24fac5e0442eef5904c09'


def replace(text, before, after):
    assert text.count(before) == 1, before
    return text.replace(before, after, 1)


def prepared_source():
    body = authority.read_text()
    body = body.replace('native-cpu-shared-value-provider', 'native-cpu-closed-dp-provider')
    body = replace(body, 'from component_gate import component, contract_inputs, terminal, HEADERS, CORES',
        'from component_gate import component, contract_inputs, terminal, HEADERS, CORES, LOCAL_CORES, DP_CALLS')
    body = replace(body, "assert contract['shared_value_interface'] == raw_component['shared_value_interface']",
        "assert contract['shared_value_interface'] == raw_component['shared_value_interface']\n" +
        "assert contract['closed_dp_interface'] == raw_component['closed_dp_interface']\n" +
        "assert contract['shared_value_interface']['raw_core_calls'] == dict(transport=2701, endpoint=3959)\n" +
        "assert contract['closed_dp_interface']['raw_core_calls'] == dict(dopri_phase=1576)")
    body = replace(body, "assert contract['basis_source_identity'] == old['full_candidate_source_identity']",
        "assert contract['basis_source_identity'] == doc(R / 'out/native-cpu-shared-value-provider/provider/source-contract.json')['full_candidate_source_identity']")
    body = replace(body, "'Sirius-native-cpu-shared-value-v45'", "'Sirius-native-cpu-closed-dp-v46'")
    body = replace(body, "'SIR-CPU-25.2.8AU'", "'SIR-CPU-25.2.8AV'")
    body = replace(body, 'assert len(parts) == 11', 'assert len(parts) == 12')
    body = replace(body, "assert len(by_name) == 11 and set(by_name) == {\n    'adapter-only.o', 'transport.o', 'endpoint.o', 'dopri_phase.o', 'optional.o', 'dispatch.o', 'v3.o',\n    'shared-dp.o', 'shared-transport.o', 'shared-endpoint.cleanup.o', 'shared-math.o'}",
        "assert len(by_name) == 12 and set(by_name) == {\n" +
        "    'adapter-only.o', 'transport.o', 'endpoint.o', 'dopri_phase.o', 'optional.o', 'dispatch.o', 'v3.o',\n" +
        "    'shared-transport.o', 'shared-endpoint.cleanup.o', 'shared-math.o', 'dp.closed.o', 'math.o'}")
    before = "for item, scope in zip(raw_v4_objects, ('transport', 'endpoint', 'dp')):\n    assert any(row['symbol'] == 'sirius_scope_' + scope + '_packet' for row in item['functions'])"
    after = """closed_dp_math_object = associated(raw_component['closed_dp_math_object'])
local_math_object = Elf(R / closed_dp_math_object['object']['path'])
assert {row['name'] for row in local_math_object.symbols if row['name'].startswith('sirius_local_value_') and row['section']} == LOCAL_CORES
for name in LOCAL_CORES:
    symbol = named(local_math_object, name)
    assert (symbol['info'] & 15) == 2 and symbol['info'] >> 4 == 1 and symbol['other'] & 3 == 2 and symbol['size']
    final = named(linked, name)
    assert final['other'] & 3 == 2 or final['info'] >> 4 == 0
    assert any(row['symbol'] == name and row['linked_address'] == final['value'] and row['bytes'] == final['size']
               for row in closed_dp_math_object['functions'])
for item, target in zip(raw_v4_objects, ('sirius_scope_transport_packet', 'sirius_scope_endpoint_packet', 'sirius_closed_dp_packet')):
    assert any(row['symbol'] == target for row in item['functions'])
assert len(raw_v4_objects[2]['functions']) == 1 and raw_v4_objects[2]['functions'][0]['symbol'] == 'sirius_closed_dp_packet'
def private_target(slug):
    return 'sirius_closed_dp_packet' if slug == 'dopri_phase' else 'sirius_scope_' + slug + '_packet'
"""
    body = replace(body, before, after)
    body = replace(body, 'for obj in raw_v4_objects + [shared_value_math_object] for row in obj[\'functions\']}',
        'for obj in raw_v4_objects + [shared_value_math_object, closed_dp_math_object] for row in obj[\'functions\']}')
    old_target = "'sirius_scope_' + ('dp' if slug == 'dopri_phase' else slug) + '_packet'"
    assert body.count(old_target) == 2
    body = body.replace(old_target, 'private_target(slug)')
    before = """for item, stage, count in zip(raw_v4_objects, ('transport', 'endpoint', 'dopri_phase'), (2701, 3959, 1581)):
    joins = direct_joins(R / item['object']['path'], item['functions'], lambda n: n.startswith('sirius_shared_value_'))
    expected = CORES - ({'sirius_shared_value_tagged_sqrt8'} if stage == 'dopri_phase' else set())
    assert len(joins) == count and {row['target'] for row in joins} == expected
    raw_direct_helper_joins[item['object']['path']] = joins
assert sum(map(len, raw_direct_helper_joins.values())) == 8241"""
    after = """for item, stage, count in zip(raw_v4_objects, ('transport', 'endpoint', 'dopri_phase'), (2701, 3959, 1576)):
    prefix = 'sirius_local_value_' if stage == 'dopri_phase' else 'sirius_shared_value_'
    joins = direct_joins(R / item['object']['path'], item['functions'], lambda n: n.startswith(prefix))
    expected = set(DP_CALLS) if stage == 'dopri_phase' else CORES
    assert len(joins) == count and {row['target'] for row in joins} == expected
    if stage == 'dopri_phase':
        from collections import Counter
        assert dict(Counter(row['target'] for row in joins)) == DP_CALLS
    raw_direct_helper_joins[item['object']['path']] = joins
assert sum(map(len, raw_direct_helper_joins.values())) == 8236"""
    body = replace(body, before, after)
    body = replace(body, "n.startswith(('sirius_simd_', 'sirius_scope_'))", "n.startswith(('sirius_simd_', 'sirius_scope_', 'sirius_closed_'))")
    body = replace(body, "shared_value_interface = dict(raw_component['shared_value_interface'], no_FP_observers=True)",
        "shared_value_interface = dict(raw_component['shared_value_interface'], no_FP_observers=True)\n" +
        "closed_dp_interface = dict(raw_component['closed_dp_interface'])\n" +
        "assert closed_dp_interface['passed'] and closed_dp_interface['no_ordinary_writes'] and closed_dp_interface['no_FP_observers']")
    body = replace(body, 'shared_value_math_object=shared_value_math_object, shared_value_interface=shared_value_interface,',
        'shared_value_math_object=shared_value_math_object, shared_value_interface=shared_value_interface,\n' +
        '    closed_dp_math_object=closed_dp_math_object, closed_dp_interface=closed_dp_interface,\n' +
        '    private_value_math_objects=[shared_value_math_object, closed_dp_math_object],')
    body = replace(body, 'eleven_input_merged_code=merge_association, eleven_input_final_code=part_association,',
        'twelve_input_merged_code=merge_association, twelve_input_final_code=part_association,')
    start = body.index("    scope='Static AU eleven-input composition only.")
    end = body.index("')\nwith (P / 'linked-native-proof.json')", start) + 2
    body = body[:start] + "    scope='Static AV twelve-part composition only. First7, both Mesa headers/cs_exec_fn, public AS/GCC/dispatcher objects and AU admitted T/E/shared9 remain exact. Only admitted DP v4 targets the one hidden six-argument closed DP entry and separate strict local-effects cc92 math9. All2701/3959 shared T/E and1576 local DP PLT32/-4 fields join actual unique hidden definitions, with18 distinct typed core definitions. Source/context confinement and complete genuine-job FP scope are bound, but actual current job PCs/results remain the required next original C8/all32 gate. No numerical transfer, old-hardware execution, ABBA/frame/scientific/release acceptance is established; ABBA/frame are deferred pending remaining-stage closure.')" + body[end:]
    return body


if __name__ == '__main__':
    exec(compile(prepared_source(), str(Path(__file__)), 'exec'))
