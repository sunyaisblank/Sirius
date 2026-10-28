"""Reuse every pinned AV code/FP/join check and prove the final AW FDE."""
from pathlib import Path
import hashlib
import importlib.util

R = Path.cwd()
authority = R / 'out/native-cpu-closed-dp-provider/verify_provider.py'
assert hashlib.sha256(authority.read_bytes()).hexdigest() == 'f88c0ad3d3584244d8b27ac632d6d41eba93361a2ba45453fa4ca4dcd8d51ff0'
spec = importlib.util.spec_from_file_location('AW_pinned_AV_verifier', authority)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def replace(text, before, after):
    assert text.count(before) == 1, before
    return text.replace(before, after, 1)


def prepared_source():
    body = module.prepared_source()
    body = body.replace('native-cpu-closed-dp-provider', 'native-cpu-closed-dp-unwind-provider')
    body = replace(body,
        "assert contract['basis_source_identity'] == doc(R / 'out/native-cpu-shared-value-provider/provider/source-contract.json')['full_candidate_source_identity']",
        "assert contract['basis_source_identity'] == doc(R / 'out/native-cpu-closed-dp-provider/provider/source-contract.json')['full_candidate_source_identity']")
    body = replace(body, "'Sirius-native-cpu-closed-dp-v46'", "'Sirius-native-cpu-closed-dp-unwind-v47'")
    body = replace(body, "'SIR-CPU-25.2.8AV'", "'SIR-CPU-25.2.8AW'")
    body = replace(body, "'dp.closed.o'", "'dp.unwind.o'")
    body = replace(body, "assert contract['closed_dp_interface'] == raw_component['closed_dp_interface']",
        "assert contract['closed_dp_interface'] == raw_component['closed_dp_interface']\n" +
        "assert contract['dp_unwind_repair'] == raw_component['dp_unwind_repair']")
    anchor = '# Count the exact retained functions in these input objects, rather than\n'
    addition = '''# Reuse the exact owned read-only FDE checker on the final DSO as well.
# Its ELF reader is replaced only by the existing byte-for-byte cached reader,
# avoiding repeated whole-section copies while preserving every check.
import importlib.util
unwind_path = R / 'out/native-cpu-closed-dp-unwind/verify.py'
unwind_plan = doc(R / 'out/native-cpu-closed-dp-unwind/build-plan.json')
checked(unwind_path, unwind_plan['inputs'][str(unwind_path.relative_to(R))])
unwind_spec = importlib.util.spec_from_file_location('AW_final_FDE_checker', unwind_path)
unwind_module = importlib.util.module_from_spec(unwind_spec)
unwind_spec.loader.exec_module(unwind_module)
unwind_module.Elf = Elf
closed_dp_unwind_table = unwind_module.verify_unwind(library, 'sirius_closed_dp_packet')
object_table = raw_component['dp_unwind_repair']['object_unwind_table']
for key in ('relative_range', 'relative_CFI_rows', 'initial_CFA', 'body_CFA',
            'saved_GP_rules', 'saved_return_PC', 'covered_call_sites', 'checked_return_offsets'):
    assert closed_dp_unwind_table[key] == object_table[key], key
for key in ('full_function_FDE', 'complete_GP_save_rules', 'shrinkwrap_and_epilogues_checked'):
    assert closed_dp_unwind_table[key] is True and object_table[key] is True, key
assert object_table['asynchronous'] is True
assert closed_dp_unwind_table['function_bytes'] == raw_v4_objects[2]['functions'][0]['bytes']
closed_dp_unwind_table.update(asynchronous=True, exact_relative_CFI_matches_repaired_object=True,
    exact_unwind_verifier=dict(path=str(unwind_path.relative_to(R)), **pin(unwind_path)),
    attribute_only_code_proof=raw_component['dp_unwind_repair']['code_proof'])
'''
    body = replace(body, anchor, addition + anchor)
    body = replace(body,
        'closed_dp_math_object=closed_dp_math_object, closed_dp_interface=closed_dp_interface,',
        "closed_dp_math_object=closed_dp_math_object, closed_dp_interface=closed_dp_interface,\n" +
        "    dp_unwind_repair=raw_component['dp_unwind_repair'], closed_dp_unwind_table=closed_dp_unwind_table,")
    body = replace(body, 'Static AV twelve-part composition only.',
        'Static AW twelve-part composition and final full-function unwind proof only. Eleven AV constituents are reused; the one repaired DP input has identical executable bytes/code relocations and new asynchronous FDE metadata. The old C8 passive observer stop is preserved and supplies no numerical verdict.')
    return body


if __name__ == '__main__':
    exec(compile(prepared_source(), str(Path(__file__)), 'exec'))
