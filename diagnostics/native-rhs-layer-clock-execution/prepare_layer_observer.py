"""Prepare an unadopted native-default Transport layer observer, without execution."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, subprocess
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
PRIOR = ROOT / 'attestations/native-vulkan/2053ef7/native-shader-clock-admission/files'
BASELINE = '2053ef79c90752d37b0075b8137cfb9b571b6820'
def seal(p):
    return {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
def one(text, old, new):
    assert text.count(old) == 1, (old, text.count(old))
    return text.replace(old, new)
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == BASELINE
assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
source = PRIOR / 'phase/retained_transport.slang'
assert seal(source) == {'bytes': 28770, 'sha256': 'b0deb2b5744899bb7e3bce8aeb81d27d0b38c81d5460ded064d48fbc1220c74e'}
compiler_receipt = json.loads((PRIOR / 'phase-compiler.json').read_text())
assert compiler_receipt['accepted']
arguments = compiler_receipt['reused_transport_compiler_run']['arguments']
assert '-DSIRIUS_RETAINED_LAYERS=117' in arguments
assert '-DSIRIUS_RETAINED_SCHWARZSCHILD_LAYERS=42' in arguments
row_words = 36 + 7 * 117 * 10
assert row_words == 8226
helper = '''
// Three lane-zero source cuts per bounded stage/layer. The original group
// join remains authoritative; these inclusive observations do not isolate
// exclusive arithmetic, physical memory stalls or ordinary benefit.
// Per-entry words: start, end, prefix-validity mask, reserved zero, three uint2 ticks.
void ObserveLayer(uint row, uint lane, uint stage, uint layer, uint slot,
                  uint start, uint end) {
    if (lane != 0) return;
    uint2 tick = getRealtimeClock();
    uint base = row * 8226u + 36u + (stage * 117u + layer) * 10u;
    phaseObservations[base + 4u + 2u * slot] = tick.x;
    phaseObservations[base + 5u + 2u * slot] = tick.y;
    if (slot == 0) {
        phaseObservations[base] = start;
        phaseObservations[base + 1u] = end;
    }
    phaseObservations[base + 2u] |= 1u << slot;
}
'''
before = source.read_text()
assert before.count('row * 36u') == 2
text = before.replace('row * 36u', 'row * 8226u')
text = one(text, '\nRetainedTriple ReadStored(uint base) {', helper + '\nRetainedTriple ReadStored(uint base) {')
pre = '''        if (start >= needed) break;
        if (start != expected'''
text = one(text, pre, '''        if (start >= needed) break;
        ObserveLayer(row, lane, stage, layer, 0, start, end);
        if (start != expected''')
post = '''#endif
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
        GroupMemoryBarrierWithGroupSync();
#endif
        if (retainedStatus != 0) return false;
        expected = end;'''
text = one(text, post, '''#endif
        ObserveLayer(row, lane, stage, layer, 1, start, end);
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
        GroupMemoryBarrierWithGroupSync();
#endif
        ObserveLayer(row, lane, stage, layer, 2, start, end);
        if (retainedStatus != 0) return false;
        expected = end;''')
# Recover the byte-exact previously reviewed observer after deleting only
# this diagnostic's helper/calls and reversing its private sidecar stride.
recovered = text.replace(helper, '')
for slot in range(3):
    call = f'        ObserveLayer(row, lane, stage, layer, {slot}, start, end);\n'
    assert recovered.count(call) == 1
    recovered = recovered.replace(call, '')
recovered = recovered.replace('row * 8226u', 'row * 36u')
assert recovered == before
phase = WORK / 'phase'
phase.mkdir()
output = phase / 'retained_transport.slang'
output.write_text(text)
inputs = [Path(__file__), source, PRIOR / 'phase-compiler.json', ROOT / 'scripts/build-retained-kernels.py']
inputs.extend((ROOT / 'src/sirius/kernels').glob('retained_*.slang'))
record = {'baseline_revision': BASELINE, 'input_seals': {str(p): seal(p) for p in inputs},
          'output': {str(output): seal(output)}, 'row_words': row_words, 'stages': 7,
          'maximum_layers': 117, 'layer_entry_words': 10, 'layer_clock_slots': 3,
          'prior_phase_cuts': 16, 'recovered_prior_observer_byte_exact': True,
          'logical_sidecar_bytes_for_two_slots': {str(cap): cap * row_words * 4 * 2 for cap in (24, 32)},
          'scope': 'Source preparation only. No compiler/device/pipeline/dispatch. Future host must admit real allocation plus original factory under unchanged8MiB, clear whole sidecars, preserve original tests and capture exact actual input tables for opcode/range joins. Inclusive perturbed lane0 cuts; no exclusive cost/ordinary benefit/scientific or release verdict.'}
(WORK / 'preparation.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'output': record['output'], 'row_words': row_words,
                  'sidecar_bytes': record['logical_sidecar_bytes_for_two_slots']}))
