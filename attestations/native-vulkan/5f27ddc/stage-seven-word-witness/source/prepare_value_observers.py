"""Prepare separate word observers without changing production kernels.

Copy only existing results at their writes. Removing the declared diagnostic
edits must recover the original source exactly. No compile/device execution.
"""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib
import json
import subprocess

ROOT = Path('/home/astra/.project/Sirius')
WORK = Path(__file__).resolve().parent
HEAD = '9fbd0e8b0d549a13476d066f2877dcd2cf9a94d9'
MAPPING = ROOT / 'attestations/review/9fbd0e8/stage7-write-mapping/files/mapping-result.json'
MAPPING_SHA = '6649a641dba12b2295476dd1951a33b7af190cd5ea2ff4a023c10cd40f8e6544'
SOURCE_SHA = {'Transport': 'd7fc2b531cf62ea3fc8454a2d1578c64c11614ac9ba046198879de925ef9f632',
              'Endpoint': '0c98ba39f13678ae4b92bd0803581379a2ded6a91bf4a6f456ee0f0cd8b5377e'}

def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

assert Path.cwd().resolve() == ROOT
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == HEAD
assert not subprocess.check_output(['git', 'status', '--porcelain'])
assert seal(MAPPING)['sha256'] == MAPPING_SHA
mapping = json.loads(MAPPING.read_text())
assert mapping['accepted_source_mapping'] and mapping['source_revision'] == HEAD
assert mapping['runtime_equality_proved'] is False
schema = {'version': 1, 'row_words': 768, 'key_base': 16, 'key_record_words': 7,
          'key_records': [0, 1, 2, 3, 4, 5, 6, 7, 44], 'arithmetic_base': 96,
          'arithmetic_record_words': 9,
          'key_kind_bits': {'normalized_admission_or_endpoint_pack': 1, 'stage6_stored_coordinate': 2, 'actual_mapped_op1_result': 4},
          'arithmetic_fields': ['presence', 'ordinal', 'destination', 'kind_shift8_or_opcode', 'hi', 'lo', 'tail', 'error', 'valid'],
          'key_fields': ['kind_presence_bits', 'input_record', 'hi', 'lo', 'tail', 'error', 'valid'],
          'header_fields': {'0': 'version', '1': 'stage1Transport2Endpoint', '2': 'offset_family1General2Schwarzschild',
                            '3': 'actual_programBase', '4': 'stored_rows', '5': 'last_successful_root_count',
                            '6': 'publication_bits1TransportStage7_2EndpointPrefix_4EndpointSuffix'},
          'scope': 'Capture only. Family derives from selected offset; canonical-table classification is deferred to complete actual input snapshots. Invalid/unvisited/flat/refused rows remain distinguishable. No result substitution.'}
schema['families'] = {}
outputs = {}
inputs = {str(MAPPING.relative_to(ROOT)): seal(MAPPING), str(Path(__file__).relative_to(ROOT)): seal(Path(__file__))}

for stage in ['Transport', 'Endpoint']:
    source = ROOT / 'src/sirius/kernels' / ('retained_' + stage.lower() + '.slang')
    assert seal(source)['sha256'] == SOURCE_SHA[stage]
    original = source.read_text()
    edits = []
    def change(old, new, count=1):
        global text
        assert text.count(old) == count, (stage, old, text.count(old), count)
        text = text.replace(old, new)
        edits.append({'old': old, 'new': new, 'count': count})
    text = original
    stage_id, stride = (1, 230) if stage == 'Transport' else (2, 225)
    families = {}
    for family, result in mapping['families'].items():
        arithmetic = [record for record in result['matches'] if record['opcode'] not in (0, 1)]
        leaves = [record for record in result['matches'] if record['opcode'] == 1]
        assert len(arithmetic) in (71, 26)
        entries = [{'slot': slot, 'ordinal': record[stage]['instruction_ordinal'],
                    'source_node': record[stage]['source_node'], 'destination': record[stage]['destination'],
                    'opcode': record['opcode'], 'inputs': record['input_records']}
                   for slot, record in enumerate(arithmetic)]
        leaf_entries = [{'ordinal': record[stage]['instruction_ordinal'], 'input': record['input_records'][0]}
                        for record in leaves]
        assert len({record['input'] for record in leaf_entries}) == len(leaf_entries)
        families[family] = {'arithmetic': entries, 'leaves': leaf_entries,
                            'selected_table_offset_words': result['programs'][stage]['selected_family_table_offset_words']}
    schema['families'][stage] = families
    def lookup(kind):
        result = []
        for family, suffix in [('general', 'General'), ('Schwarzschild', 'Schwarzschild')]:
            cases = families[family][kind]
            result += [f'uint WordLookup{kind.title()}{suffix}(uint index) {{', '    switch (index) {']
            for record in cases:
                result.append(f'    case {record["ordinal"]}u: return {record["slot"] if kind == "arithmetic" else record["input"]}u;')
            result += ['    default: return 0xffffffffu;', '    }', '}']
        return '\n'.join(result)
    helper = '''
// BEGIN PRIVATE WORD OBSERVER: no original arithmetic, writes or joins removed.
[[vk::binding(2, 0)]] RWStructuredBuffer<uint> wordObservations;
bool WordRowAvailable(uint row) {
    uint words, stride;
    wordObservations.GetDimensions(words, stride);
    return row < words / 768u;
}
uint WordFamily(uint programBase) {
    uint originalBase = 1u + STRIDEu * inputs[0];
    if (programBase == originalBase) return 1u;
    if (programBase == originalBase + SIRIUS_RETAINED_GENERAL_PROGRAM_WORDS) return 2u;
    return 0u;
}
void WordCopy(uint base, RetainedTriple value) {
    wordObservations[base] = RPBits(value.hi);
    wordObservations[base + 1u] = RPBits(value.lo);
    wordObservations[base + 2u] = RPBits(value.tail);
    wordObservations[base + 3u] = RPBits(value.error);
    wordObservations[base + 4u] = value.valid;
}
void WordKey(uint row, uint inputIndex, RetainedTriple value, uint kind) {
    if (!WordRowAvailable(row) || (inputIndex > 7u && inputIndex != 44u)) return;
    uint slot = inputIndex == 44u ? 8u : inputIndex;
    uint base = row * 768u + 16u + slot * 7u;
    wordObservations[base + 1u] = inputIndex;
    WordCopy(base + 2u, value);
    uint previous;
    InterlockedOr(wordObservations[base], kind, previous);
}
LOOKUPS
void WordEvaluated(uint row, uint programBase, uint index, uint destination,
                   uint op, uint a, RetainedTriple value, uint kind) {
    if (!WordRowAvailable(row)) return;
    TRANSPORT_STAGE_GUARD
    uint family = WordFamily(programBase);
    if (family == 0u) return;
    uint leaf = family == 1u ? WordLookupLeavesGeneral(index) : WordLookupLeavesSchwarzschild(index);
    // Only the original unique mapped leaf write may update its key. Altered
    // tables cannot turn arbitrary simultaneous writes into key writers.
    if (leaf != 0xffffffffu && op == 1u && a == leaf) WordKey(row, leaf, value, 4u);
    uint slot = family == 1u ? WordLookupArithmeticGeneral(index) : WordLookupArithmeticSchwarzschild(index);
    if (slot == 0xffffffffu) return;
    uint base = row * 768u + 96u + slot * 9u;
    wordObservations[base + 1u] = index;
    wordObservations[base + 2u] = destination;
    wordObservations[base + 3u] = (kind << 8u) | op;
    WordCopy(base + 4u, value);
    wordObservations[base] = 1u; // Presence is separate from numerical valid0.
}
void WordPublished(uint row, uint programBase, uint roots, uint bits) {
    if (!WordRowAvailable(row)) return;
    uint base = row * 768u;
    wordObservations[base] = 1u;
    wordObservations[base + 1u] = STAGEu;
    wordObservations[base + 2u] = WordFamily(programBase);
    wordObservations[base + 3u] = programBase;
    wordObservations[base + 4u] = inputs[0];
    wordObservations[base + 5u] = roots;
    uint previous;
    InterlockedOr(wordObservations[base + 6u], bits, previous);
}
// END PRIVATE WORD OBSERVER.
'''.replace('STRIDEu', str(stride) + 'u').replace('STAGEu', str(stage_id) + 'u')
    helper = helper.replace('LOOKUPS', lookup('arithmetic') + '\n' + lookup('leaves'))
    helper = helper.replace('TRANSPORT_STAGE_GUARD', 'if (kind != 6u) return;' if stage == 'Transport' else '')
    # Replace STAGE only in its literal metadata placeholder, not the guard name.
    assert 'TRANSPORT_' not in helper
    change('\nRetainedTriple ReadStored(uint base) {', helper + '\nRetainedTriple ReadStored(uint base) {')
    change('uint programBase, uint registers) {', 'uint programBase, uint registers, uint wordKind) {')
    change('    WriteValue(destination, r, rowBase);\n',
           '    WriteValue(destination, r, rowBase);\n    WordEvaluated(row, programBase, index, destination, op, a, r, wordKind);\n')
    kind = 'stage' if stage == 'Transport' else 'required'
    change('EvaluateOne(node, row, rowBase, programBase, registers)', f'EvaluateOne(node, row, rowBase, programBase, registers, {kind})', 2)
    change('EvaluateOne(start + lane, row, rowBase, programBase, registers)', f'EvaluateOne(start + lane, row, rowBase, programBase, registers, {kind})')
    if stage == 'Transport':
        change('    for (uint i = 0; i < 46; ++i)\n        if (ReadOriginal(i, row).valid == 0) return;',
               '''    for (uint i = 0; i < 46; ++i) {
        RetainedTriple wordOriginal = ReadOriginal(i, row);
        if (lane == 0 && (i < 4 || i == 44)) WordKey(row, i, wordOriginal, 1u);
        if (wordOriginal.valid == 0) return;
    }''')
        change('            WriteStored(rowBase + 804 + i * 5, value);',
               '            WriteStored(rowBase + 804 + i * 5, value);\n            if (stage == 6 && i < 4) WordKey(row, 4 + i, value, 2u);')
        change('    return true;\n}\nRetainedTriple Rational',
               '    if (lane == 0 && stage == 6) WordPublished(row, programBase, 40u, 1u);\n    return true;\n}\nRetainedTriple Rational')
    else:
        change('        if (lane == 0) WriteStored(rowBase + 504 + 5 * i, value);',
               '        if (lane == 0) WriteStored(rowBase + 504 + 5 * i, value);\n        if (lane == 0) WordKey(row, i, value, 1u);', 2)
        change('    return true;\n}\n\nRetainedTriple Metric',
               '    if (lane == 0) WordPublished(row, programBase, required, required == 20 ? 2u : 4u);\n    return true;\n}\n\nRetainedTriple Metric')
    recovered = text
    for edit in reversed(edits):
        assert recovered.count(edit['new']) == edit['count'], (stage, edit['old'])
        recovered = recovered.replace(edit['new'], edit['old'])
    assert recovered == original
    output = WORK / 'kernels' / source.name
    output.parent.mkdir(parents=True, exist_ok=True)
    assert not output.exists()
    output.write_text(text)
    outputs[str(output.relative_to(ROOT))] = seal(output)
    inputs[str(source.relative_to(ROOT))] = seal(source)
    schema['families'][stage]['source_edits'] = edits

for path in (ROOT / 'src/sirius/kernels').glob('retained_*.slang'):
    inputs[str(path.relative_to(ROOT))] = seal(path)
for name, value in [('value-schema.json', schema), ('preparation.json', {
    'accepted_source_preparation': True, 'source_revision': HEAD, 'input_seals': inputs, 'outputs': outputs,
    'production_sources_recovered_byte_exact': True, 'row_words': 768,
    'logical_two_sidecar_bytes': {str(capacity): capacity * 768 * 4 * 2 for capacity in [24, 32]},
    'scope': 'Private source-only diagnostic preparation. No compiler/device/candidate. Original calls/writes/guards/barriers retained; normalization calls materialized once without extra arithmetic. Full actual T/E input/output plus per-command cleared sidecars required for execution.'})]:
    path = WORK / name
    assert not path.exists()
    path.write_text(json.dumps(value, indent=2) + '\n')
print(json.dumps({'outputs': outputs, 'row_words': 768, 'source_preparation_only': True}))
