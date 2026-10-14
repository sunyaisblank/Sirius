"""Root-operated source-reader controls on two complete immutable shader dumps."""
from pathlib import Path
import hashlib, json
from multiply_profile import function_body, multiply_shape

W = Path(__file__).resolve().parent
P = json.loads((W / 'packet.json').read_text())

def fixture(name):
    item = P['parser_fixtures'][name]
    data = Path(item['path']).read_bytes()
    assert len(data) == item['bytes'] and hashlib.sha256(data).hexdigest() == item['sha256']
    return data.decode()

native, ssa = fixture('native_o0'), fixture('portable_ssa')
body = function_body(native, 'PB32Multiply')

def once(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new, 1)

def changed(old, new):
    updated = once(body, old, new)
    return once(native, body, updated)

swapped = once(body, '               OpStore %shift %uint_21\n', '               OpStore %shift %uint_20\n')
# Disambiguate by the exact false arm; the true arm was already changed.
swapped = once(swapped, '%1806 = OpLabel\n               OpStore %shift %uint_20',
               '%1806 = OpLabel\n               OpStore %shift %uint_21')
rows = [
    ('actual_complete_native_o0', native, True, True),
    ('accepted_complete_portable_ssa', ssa, False, True),
    ('o0_requires_explicit_admission', native, False, False),
    ('reversed_21_20_arms', once(native, body, swapped), True, False),
    ('wrong_selector_predicate', changed('OpBranchConditional %upper %1805 %1806',
                                       'OpBranchConditional %1795 %1805 %1806'), True, False),
    ('wrong_true_constant', changed('               OpStore %shift %uint_21\n', '               OpStore %shift %uint_1\n'), True, False),
    ('intervening_distance_store', changed('%1810 = OpLoad %uint %shift',
        '%o0_store_marker = OpCopyObject %uint %uint_20\n               OpStore %shift %o0_store_marker\n       %1810 = OpLoad %uint %shift'), True, False),
    ('pointer_copy_alias', changed('%1760 = OpVariable %_ptr_Function_int Function',
        '%1760 = OpVariable %_ptr_Function_int Function\n       %o0_alias = OpCopyObject %_ptr_Function_uint %shift'), True, False),
    ('pointer_call_escape', changed('%1760 = OpVariable %_ptr_Function_int Function',
        '%1760 = OpVariable %_ptr_Function_int Function\n       %o0_escape = OpFunctionCall %bool %PB32Zero %shift'), True, False),
    ('true_arm_side_entrance', changed('OpBranchConditional %1795 %1797 %1796',
                                    'OpBranchConditional %1795 %1805 %1796'), True, False),
    ('mask_distance_mismatch', changed('OpShiftLeftLogical %ulong %ulong_1 %1807',
                                     'OpShiftLeftLogical %ulong %ulong_1 %uint_1'), True, False),
    ('mask_load_from_sticky_epoch', changed('OpShiftLeftLogical %ulong %ulong_1 %1807',
                                          'OpShiftLeftLogical %ulong %ulong_1 %1817'), True, False),
    ('sticky_uses_distance', changed('OpBitwiseOr %uint %1812 %1817',
                                   'OpBitwiseOr %uint %1812 %1810'), True, False),
    ('sticky_wrong_predicate', changed('OpBranchConditional %1813 %1815 %1816',
                                      'OpBranchConditional %upper %1815 %1816'), True, False),
    ('missing_true_store', changed('               OpStore %shift %uint_21\n', ''), True, False),
]
results = []
for name, text, allow_o0, expected in rows:
    try:
        proof = multiply_shape(text, allow_o0=allow_o0)
        actual, error = True, None
    except (AssertionError, ValueError, KeyError, IndexError) as failure:
        actual, error = False, type(failure).__name__
    assert actual == expected, (name, expected, actual)
    if name == 'actual_complete_native_o0':
        choices = proof['function_slot_load_choices']
        assert len(choices) == 3 and len({item['merge'] for item in choices}) == 2
        assert {item['load'] for item in choices} == {'%1807', '%1810', '%1817'}
    if name == 'accepted_complete_portable_ssa':
        assert proof['function_slot_load_choices'] == []
    results.append({'name': name, 'expected_acceptance': expected, 'actual_acceptance': actual,
                    'failure_type': error, 'full_fixture_sha256': hashlib.sha256(text.encode()).hexdigest()})
target = W / 'parser-controls-result.json'
assert not target.exists()
target.write_text(json.dumps({'passed': True, 'controls': results,
    'scope': 'Full immutable dump source-reader controls only. Some negative textual perturbations intentionally violate SPIR-V SSA/type/CFG validity; none is compiled, validated, dispatched, or a numerical fixture.'},
    sort_keys=True, indent=2) + '\n')
print('Complete-module reader controls passed; no shader/tool/provider executed.', flush=True)
