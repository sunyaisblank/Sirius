"""Retain the pinned AO ELF/ISA proof; bind only the two AP routing objects."""
from pathlib import Path
from authority import source, replace
from component_gate import component, contract_inputs, doc, pin, terminal

new_objects = component()
P = Path(__file__).resolve().parent / 'provider'
contract = doc(P / 'source-contract.json')
contract_inputs(contract)
terminal(P / 'build/owner.json')
assert contract['guard_shared_byte_authority'] == 'out/native-cpu-endpoint-header-domain/native_packet.hpp'
assert contract['source_inputs'][contract['guard_shared_byte_authority']] == pin(Path.cwd() / contract['guard_shared_byte_authority'])
old = doc(Path.cwd() / 'out/native-cpu-v4-provider/provider/source-contract.json')
assert contract['adapter_input_objects'] == old['adapter_input_objects']
assert contract['component_manifest'] == old['component_manifest']
expected = {Path(row['path']).name: row for row in new_objects}
assert contract['native_objects'][1] == expected['endpoint.o']
assert contract['native_objects'][3] == expected['admission.o']

body = source('out/native-cpu-v4-provider/verify_provider.py', 'e60ea2063152bc100a5e738bf73fa335d8f44380385f73fd67f7d24a54b3335b')
body = replace(body, "assert W == R / 'out/native-cpu-v4-provider'", "assert W == R / 'out/native-cpu-endpoint-header-provider'")
body = replace(body, "old = doc(R / 'out/native-cpu-local-identities/provider/source-contract.json')",
    "old = doc(R / 'out/native-cpu-v4-provider/provider/source-contract.json')")
body = replace(body, "assert len(objects) == 8 and objects[:7] == old['native_objects'][:7]",
    "assert len(objects) == 8 and all(objects[i] == old['native_objects'][i] for i in (0, 2, 4, 5, 6, 7))")
body = replace(body, "contract['new_native_compilations'] == 0", "contract['new_native_compilations'] == 2")
body = replace(body, "    first7_native_objects_unchanged=True, original_v3_component_proof_reused=True,",
    "    unchanged_AO_native_object_indices=[0, 2, 4, 5, 6, 7], changed_AO_native_object_indices=[1, 3],\n"
    "    only_Endpoint_scalar_and_pure_admission_replaced=True, original_v3_component_proof_reused=True,\n"
    "    guard_shared_byte_authority=contract['guard_shared_byte_authority'],\n"
    "    header_domain_compatibility=pin(R / contract['header_domain_component']['compatibility']),")
body = replace(body, "with (P / 'linked-native-proof.json').open('x') as stream:",
    "assert sum(len(row['functions']) for row in rows) == 920\n"
    "with (P / 'linked-native-proof.json').open('x') as stream:")
exec(compile(body, str(Path(__file__)), 'exec'))
