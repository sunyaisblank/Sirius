"""Prepare exact-order native division observations; no compiler or device work."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, re, struct, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
BASELINE = '2053ef79c90752d37b0075b8137cfb9b571b6820'
PRIOR = ROOT / 'attestations/native-vulkan/2053ef7/native-shader-clock-admission/files'
LAYERS = ROOT / 'attestations/native-vulkan/24a8e14/native-rhs-layer-clock-execution/files'
STRIDE = 8226
SLOTS = 12
ENTRY_WORDS = 56

def read(path):
    return json.loads(path.read_text())

def seal(path):
    return {'bytes': path.stat().st_size,
            'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}

def git(*args):
    return subprocess.check_output(['git', *args], text=True, timeout=30).strip()

def one(text, old, new):
    assert text.count(old) == 1, (old, text.count(old))
    return text.replace(old, new)

assert WORK == ROOT / 'out/native-division-phase-observer'
assert git('rev-parse', 'HEAD') == BASELINE == git('rev-parse', '@{u}')
assert not git('status', '--porcelain')
original = PRIOR / 'phase/retained_transport.slang'
assert seal(original) == {'bytes': 28770, 'sha256': 'b0deb2b5744899bb7e3bce8aeb81d27d0b38c81d5460ded064d48fbc1220c74e'}
analysis_path = LAYERS / 'device-timestamps/observation-analysis.json'
summary = read(LAYERS / 'results.json')
assert summary['accepted'] and seal(analysis_path) == summary['cases']['device-timestamps']['analysis']
analysis = read(analysis_path)
assert analysis['accepted'] and analysis['source_revision'] == '24a8e145d1da5dca75a4fa26ec8e639020eecf37'
header_path = ROOT / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
assert seal(header_path) == {'bytes': 7952166, 'sha256': '3f7aca3372bb872934e978fb0270f2375ce8d2361555b35c8b891e42e07161cc'}
program_matches = re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> kTransportProgram\{\{(.*?)\}\};', header_path.read_text(), re.S)
assert len(program_matches) == 1
count, body = program_matches[0]
words = tuple(int(x[:-1], 0) for x in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u', body))
assert len(words) == int(count) == 15967
selections = {}
for family, offset, length, nodes_count, layer_count in (
        ('general', 0, 13096, 2587, 117), ('schwarzschild', 13096, 2871, 557, 42)):
    tables = [t for t in analysis['actual_tables'].values() if t['family'] == family]
    assert len(tables) == 1
    table = tables[0]
    data = struct.pack('<' + 'I' * length, *words[offset:offset + length])
    assert table['key'] == family + ':' + hashlib.sha256(data).hexdigest()
    assert table['count'] == nodes_count and table['layers'] == layer_count
    selected = []
    for layer, (start, end) in enumerate(zip(table['offsets'], table['offsets'][1:])):
        if not all(node[0] == 5 for node in table['nodes'][start:end]):
            continue
        observed = [r for r in analysis['published_layer_summary']
                    if r['table'] == table['key'] and r['layer'] == layer]
        assert len(observed) == 1 and observed[0]['opcode_composition'] == {'5': end - start}
        selected.append({'slot': len(selected), 'layer': layer, 'node': start,
                         'layer_end': end, 'node_words': table['nodes'][start],
                         'prior_inclusive_body': observed[0]['lane_zero_before_to_after_body']})
    assert len(selected) == (12 if family == 'general' else 4)
    selections[family] = {'program_words': nodes_count, 'table_key': table['key'], 'nodes': selected}

selector = 'uint DivisionSlot(uint index, uint programBase) {\n    uint count = inputs[programBase];\n'
for family, selection in selections.items():
    selector += f"    if (count == {selection['program_words']}u) {{\n"
    for node in selection['nodes']:
        selector += f"        if (index == {node['node']}u) return {node['slot']}u;\n"
    selector += '    }\n'
selector += '    return 12u;\n}\n'
helpers = '''
// Private observations replace the earlier layer records while retaining its
// allocation stride. Each selected lane has a disjoint stage/node entry.
// Metadata: node, program count, reached mask, terminal reason, returned valid,
// magic; exact a/b five-word snapshots; then twenty uint2 device-clock samples.
uint DivisionBase(uint row, uint stage, uint slot) {
    return row * 8226u + 36u + (stage * 12u + slot) * 56u;
}
void StartDivision(uint row, uint stage, uint slot, uint index, uint count,
                   RetainedTriple a, RetainedTriple b) {
    if (slot >= 12u) return;
    uint base = DivisionBase(row, stage, slot);
    phaseObservations[base] = index;
    phaseObservations[base + 1u] = count;
    phaseObservations[base + 5u] = 0x4456434bu;
    phaseObservations[base + 6u] = RPBits(a.hi);
    phaseObservations[base + 7u] = RPBits(a.lo);
    phaseObservations[base + 8u] = RPBits(a.tail);
    phaseObservations[base + 9u] = RPBits(a.error);
    phaseObservations[base + 10u] = a.valid;
    phaseObservations[base + 11u] = RPBits(b.hi);
    phaseObservations[base + 12u] = RPBits(b.lo);
    phaseObservations[base + 13u] = RPBits(b.tail);
    phaseObservations[base + 14u] = RPBits(b.error);
    phaseObservations[base + 15u] = b.valid;
}
void ObserveDivision(uint row, uint stage, uint slot, uint cut) {
    if (slot >= 12u) return;
    uint2 tick = getRealtimeClock();
    uint base = DivisionBase(row, stage, slot);
    phaseObservations[base + 16u + 2u * cut] = tick.x;
    phaseObservations[base + 17u + 2u * cut] = tick.y;
    phaseObservations[base + 2u] |= 1u << cut;
}
RetainedTriple FinishDivision(uint row, uint stage, uint slot, uint reason,
                              RetainedTriple result) {
    if (slot < 12u) {
        uint base = DivisionBase(row, stage, slot);
        phaseObservations[base + 3u] = reason;
        phaseObservations[base + 4u] = result.valid;
    }
    return result;
}
'''
triple_path = ROOT / 'src/sirius/kernels/retained_triple.slang'
match = re.search(r'RetainedTriple RTDivide\(RetainedTriple a,RetainedTriple b\) \{\n#ifdef SIRIUS_RETAINED_PORTABLE\n.*?#else\n(.*?)#endif\n\}', triple_path.read_text(), re.S)
assert match
native_body = match[1]
observed_body = native_body
changes = []
def change(old, new):
    global observed_body
    observed_body = one(observed_body, old, new)
    changes.append((old, new))
change('    if(a.valid==0||b.valid==0) return RTInvalid();',
       '    StartDivision(row,stage,slot,index,count,a,b);\n    ObserveDivision(row,stage,slot,0u);\n    if(a.valid==0||b.valid==0) return FinishDivision(row,stage,slot,2u,RTInvalid());')
change('RPFromBits(0x5e000000u))) return RTInvalid();',
       'RPFromBits(0x5e000000u))) return FinishDivision(row,stage,slot,3u,RTInvalid());')
change('    if(RTExact(b,RPOne())) return a;',
       '    if(RTExact(b,RPOne())) return FinishDivision(row,stage,slot,4u,a);')
change('    RetainedTriple q=RTConstant(RPDivideSeed(a.hi,b.hi));',
       '    RetainedTriple q=RTConstant(RPDivideSeed(a.hi,b.hi));\n    ObserveDivision(row,stage,slot,1u);')
change('        RetainedTriple residual=RTSubtract(a,RTMultiply(b,q));',
       '        RetainedTriple product=RTMultiply(b,q);\n        ObserveDivision(row,stage,slot,2u+5u*i);\n        RetainedTriple residual=RTSubtract(a,product);\n        ObserveDivision(row,stage,slot,3u+5u*i);')
change('            return RTInvalid();',
       '            return FinishDivision(row,stage,slot,5u,RTInvalid());\n        ObserveDivision(row,stage,slot,4u+5u*i);')
change('        q=RTCenter(RTAdd(q,RTConstant(RPDivideSeed(residual.hi,b.hi))));',
       '        RetainedTriple correction=RTConstant(RPDivideSeed(residual.hi,b.hi));\n        ObserveDivision(row,stage,slot,5u+5u*i);\n        q=RTCenter(RTAdd(q,correction));\n        ObserveDivision(row,stage,slot,6u+5u*i);')
change('        if(q.valid==0) return q;',
       '        if(q.valid==0) return FinishDivision(row,stage,slot,6u,q);')
change('    RetainedTriple residual=RTSubtract(a,RTMultiply(b,q));',
       '    RetainedTriple product=RTMultiply(b,q);\n    ObserveDivision(row,stage,slot,17u);\n    RetainedTriple residual=RTSubtract(a,product);\n    ObserveDivision(row,stage,slot,18u);')
change('    if(residual.valid==0) return residual;',
       '    if(residual.valid==0) return FinishDivision(row,stage,slot,7u,residual);')
change('    return RTPack(q.hi,q.lo,q.tail,RPRatioUpper(RTAbsoluteUpper(residual),lower));',
       '    RetainedTriple result=RTPack(q.hi,q.lo,q.tail,RPRatioUpper(RTAbsoluteUpper(residual),lower));\n    ObserveDivision(row,stage,slot,19u);\n    return FinishDivision(row,stage,slot,result.valid==0?8u:1u,result);')
recovered = observed_body
for old, new in reversed(changes):
    recovered = one(recovered, new, old)
assert recovered == native_body
function = ('RetainedTriple ObservedRTDivide(RetainedTriple a,RetainedTriple b,uint row,\n'
            '                                uint stage,uint index,uint slot,uint count) {\n'
            + observed_body + '}\n')

before = original.read_text()
text = before.replace('row * 36u', 'row * 8226u')
assert before.count('row * 36u') == 2
text = one(text, '\nRetainedTriple ReadStored(uint base) {',
           '\n' + selector + helpers + function + '\nRetainedTriple ReadStored(uint base) {')
text = one(text, 'bool EvaluateOne(uint index, uint row, uint rowBase, uint programBase, uint registers) {',
           'bool EvaluateOne(uint index, uint row, uint rowBase, uint programBase, uint registers, uint stage) {')
text = one(text, '                r = RTDivide(av, bv);',
           '                r = ObservedRTDivide(av,bv,row,stage,index,DivisionSlot(index,programBase),inputs[programBase]);')
assert text.count('EvaluateOne(node, row, rowBase, programBase, registers)') == 2
text = text.replace('EvaluateOne(node, row, rowBase, programBase, registers)',
                    'EvaluateOne(node, row, rowBase, programBase, registers, stage)')
text = one(text, 'EvaluateOne(start + lane, row, rowBase, programBase, registers)',
           'EvaluateOne(start + lane, row, rowBase, programBase, registers, stage)')
assert 'GroupMemoryBarrier' not in function + helpers + selector
assert text.count('getRealtimeClock()') == 2
phase = WORK / 'phase'
phase.mkdir()
target = phase / 'retained_transport.slang'
target.write_text(text)
inputs = [Path(__file__), original, PRIOR / 'phase-compiler.json', analysis_path,
          LAYERS / 'results.json', header_path, triple_path,
          ROOT / 'scripts/build-retained-kernels.py']
inputs.extend((ROOT / 'src/sirius/kernels').glob('retained_*.slang'))
record = {'baseline_revision': BASELINE, 'input_seals': {str(p): seal(p) for p in inputs},
          'output': {str(target): seal(target)}, 'selections': selections,
          'native_division_body_restored_byte_exact': True,
          'row_words': STRIDE, 'selected_slots_per_stage': SLOTS, 'entry_words': ENTRY_WORDS,
          'used_row_words': 36 + 7 * SLOTS * ENTRY_WORDS,
          'cut_layout': ['entry', 'initial_guard_and_seed_complete']
                        + [f'correction_{i}_{name}' for i in range(3)
                           for name in ('product', 'subtract', 'refusal_check', 'seed_constant', 'center_add')]
                        + ['final_product', 'final_subtract', 'final_pack'],
          'terminal_reasons': {'1': 'final_pack_valid', '2': 'invalid_operand',
                               '3': 'initial_lower_or_ratio_refusal', '4': 'unit_denominator',
                               '5': 'correction_residual_or_ratio_refusal', '6': 'invalid_updated_q',
                               '7': 'invalid_final_residual', '8': 'invalid_final_pack'},
          'scope': 'Diagnostic source preparation only. First node of every prior division-only layer, all rows and seven stages; actual a/b five-word snapshots. Existing8226stride and8MiB admission unchanged. Original guards/short-circuit/three corrections/final residual/return words retained. No new collectives. Actual compiled driver placement and original numerical/stream/refusal controls remain required; clocks do not guarantee code-motion barriers. No exclusive operation/total division cost, ordinary benefit, frame or qualification claim.'}
assert record['used_row_words'] == 4740 < STRIDE and len(record['cut_layout']) == 20
(WORK / 'preparation.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'output': record['output'], 'used_row_words': record['used_row_words'],
                  'selected_nodes': {k: len(v['nodes']) for k, v in selections.items()}}))
