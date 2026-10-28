"""Finite static gate for newly compiled baseline component/diagnostic objects."""
from pathlib import Path
import hashlib
import json
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-dispatch-admission-token'
sys.dont_write_bytecode = True
sys.path.insert(0, str(R / 'out/native-cpu-field-compatibility'))
from elf_proof import Elf


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


finite = json.loads((W / 'finite-inputs.json').read_text())
authority = R / 'out/native-cpu-inline-private-bank/check_code.py'
assert pin(authority) == finite['source_inputs'][str(authority.relative_to(R))]
text = authority.read_text()
anchor = "proof=doc(W/'derivation.json');built=doc(W/'build-result.json')"
assert text.count(anchor) == 1
prefix = text[:text.index(anchor)]
production_scope = {'__file__': str(authority), '__name__': 'sealed_AS_encoding_reader'}
exec(compile(prefix, str(authority), 'exec'), production_scope)
# The diagnostic wrapper only READS MXCSR, as the original fixtures do. No
# production FP-state operation is added or allowed, and no write is exempted.
assert prefix.count("'stmxcsr',") == 1
diagnostic_scope = {'__file__': str(authority), '__name__': 'diagnostic_encoding_reader'}
exec(compile(prefix.replace("'stmxcsr',", '', 1), str(authority), 'exec'), diagnostic_scope)
rows = []
for role in ('transport', 'endpoint', 'dopri_phase', 'rows_adapter', 'fixture_routes'):
    path = B / (role + '.o')
    scope = diagnostic_scope if role == 'fixture_routes' else production_scope
    elf, instructions, sections = scope['inspect'](path)
    assert all(row['encoding_form'][1] == 'legacy' for row in instructions), role
    definitions = {s['name']: s for s in elf.symbols if s['section'] and (s['info'] & 15) == 2}
    edges = {}
    for name, symbol in definitions.items():
        edges[name] = sorted({target['name'] for section, offset, kind, addend, target in elf.relocations
                              if section == symbol['section'] and
                              symbol['value'] <= offset < symbol['value'] + symbol['size'] and target['name']})
    if role in ('transport', 'endpoint', 'dopri_phase'):
        public = 'sirius_native_' + role + '_try_packet'
        private = 'sirius_native_' + role + '_admitted_packet'
        assert public in definitions and private in definitions
        assert definitions[private]['other'] & 3 == 2
        assert any(s['name'] == 'memcmp' and not s['section'] for s in elf.symbols)
        assert 'memcmp' not in edges[private] and 'bcmp' not in edges[private]
        assert public not in edges[private]
    elif role == 'rows_adapter':
        for stage in ('transport', 'endpoint', 'dopri_phase'):
            public = 'sirius_native_' + stage + '_rows_try_packet'
            private = 'sirius_native_' + stage + '_rows_admitted_packet'
            assert public in definitions and private in definitions
            assert definitions[private]['other'] & 3 == 2
            assert 'sirius_native_' + stage + '_admit_packet' in edges[public]
            assert 'sirius_native_' + stage + '_admit_packet' not in edges[private]
            assert 'sirius_native_' + stage + '_try_packet' not in edges[private]
            assert 'sirius_native_' + stage + '_admitted_packet' in edges[private]
            assert not {'memcmp', 'bcmp'} & set(edges[private])
    rows.append(dict(role=role, object=dict(path=str(path.relative_to(R)), **pin(path)),
        executable_sections=sections, instructions=len(instructions), functions=len(definitions),
        baseline_legacy_encodings_only=True, no_FMA_x87_or_FP_control_writes=True,
        new_export_edges={name: targets for name, targets in edges.items()
                          if name.startswith(('sirius_native_', 'sirius_token_'))}))

diagnostic = Elf(B / 'fixture_routes.o')
targets = {target['name'] for _, _, _, _, target in diagnostic.relocations}
expected = {
    'sirius_native_transport_admit_packet', 'sirius_native_endpoint_admit_packet',
    'sirius_native_transport_rows_try_packet', 'sirius_native_endpoint_rows_try_packet',
    'sirius_native_transport_rows_admitted_packet', 'sirius_native_endpoint_rows_admitted_packet',
    'sirius_native_endpoint_try_packet', 'sirius_native_endpoint_admitted_packet',
}
assert expected <= targets
result = dict(passed=True, objects=rows, source_encoding_reader=dict(path=str(authority.relative_to(R)), **pin(authority)),
    public_guards_still_reference_actual_admission_or_memcmp=True,
    private_export_direct_edges_skip_program_comparisons_only=True,
    diagnostic_real_preflight_and_private_symbol_references=sorted(expected),
    limits='Exact original baseline compiler flags plus complete emitted-section coverage and rejection of VEX/EVEX/FMA/x87/FP-state writes. The only FP-state instruction exemption is diagnostic stmxcsr read. Direct-edge checks complement source inverses; complete call-target resolution is required at link. This is not execution, proof of every path, actual Mesa capability selection, performance or qualification.')
with (W / 'code-proof.json').open('x') as stream:
    json.dump(result, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(passed=True, production_objects=4, diagnostic_objects=1)))
