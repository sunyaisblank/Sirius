"""Build unadopted native-default shader copies with existing-join observations."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json
R=Path.cwd();W=Path(__file__).resolve().parent;D=W/'phase';D.mkdir()
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def one(text,old,new):
 assert text.count(old)==1,(old,text.count(old));return text.replace(old,new)
record={'input_seals':{},'outputs':{},'scope':'Unadopted native-default source copies only. Existing joins/calculations/outputs preserved. Sidecars require real descriptor binding2, clearing and allocation accounting before execution; not yet connected.'}
for stem,words,slots in [('retained_transport',36,16),('retained_endpoint',12,4)]:
 p=R/'src/sirius/kernels'/(stem+'.slang');text=p.read_text();record['input_seals'][str(p)]=seal(p)
 helper='''
#if defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_FP64) || defined(SIRIUS_RETAINED_FMA32)
#error This diagnostic observes only the original native-default binary32 route.
#endif
// Diagnostic descriptor only. Host must clear the entire sidecar per command.
// Metadata: magic, row, cut-validity bitmap, successful completion; then uint2 ticks.
[[vk::binding(2, 0)]] RWStructuredBuffer<uint> phaseObservations;
void ObservePhase(uint row, uint lane, uint slot) {
    if (lane != 0) return;
    uint2 tick = getRealtimeClock();
    uint base = row * '''+str(words)+'''u;
    phaseObservations[base + 4 + 2 * slot] = tick.x;
    phaseObservations[base + 5 + 2 * slot] = tick.y;
    phaseObservations[base] = 0x53434c4bu;
    phaseObservations[base + 1] = row;
    phaseObservations[base + 2] |= 1u << slot;
}
'''
 text=one(text,'[[vk::binding(1, 0)]] RWStructuredBuffer<uint> outputs;', '[[vk::binding(1, 0)]] RWStructuredBuffer<uint> outputs;\n'+helper)
 if stem=='retained_transport':
  text=one(text,'    return true;\n}\nRetainedTriple Rational(', '    ObservePhase(row, lane, 2 + 2 * stage);\n    return true;\n}\nRetainedTriple Rational(')
  text=one(text,'        GroupMemoryBarrierWithGroupSync();\n#endif\n        for (uint i = lane; i < 40;', '        GroupMemoryBarrierWithGroupSync();\n#endif\n        if (stage == 0) ObservePhase(row, lane, 0);\n        for (uint i = lane; i < 40;')
  text=one(text,'        if (retainedStatus != 0) return;\n        if (lane == 0) outputs[rowBase + 1] = stage + 1;', '        ObservePhase(row, lane, 1 + 2 * stage);\n        if (retainedStatus != 0) return;\n        if (lane == 0) outputs[rowBase + 1] = stage + 1;')
  text=one(text,'    if (retainedStatus != 0) return;\n    if (lane == 0) outputs[rowBase] = 1;', '    ObservePhase(row, lane, 15);\n    if (retainedStatus != 0) return;\n    if (lane == 0) outputs[rowBase] = 1;')
  text=one(text,'    if (lane == 0) outputs[rowBase + 2] = 1;\n}', '    if (lane == 0) outputs[rowBase + 2] = 1;\n    if (lane == 0) phaseObservations[row * 36u + 3] = 1;\n}')
 else:
  text=one(text,'    GroupMemoryBarrierWithGroupSync();\n#endif\n    uint layerBase = programBase + 103', '    GroupMemoryBarrierWithGroupSync();\n#endif\n    if (required == 20) ObservePhase(row, lane, 0);\n    uint layerBase = programBase + 103')
  text=one(text,'    return true;\n}\n\nRetainedTriple Metric(', '    ObservePhase(row, lane, required == 20 ? 1 : 3);\n    return true;\n}\n\nRetainedTriple Metric(')
  text=one(text,'#endif\n    if (!EvaluateProgram(row, rowBase, programBase, count, registers, 100, lane, layers, prefix)) return;', '#endif\n    ObservePhase(row, lane, 2);\n    if (!EvaluateProgram(row, rowBase, programBase, count, registers, 100, lane, layers, prefix)) return;')
  text=one(text,'    if (lane == 0) outputs[rowBase + 1] = component;\n}\n#endif', '    if (lane == 0) outputs[rowBase + 1] = component;\n    if (lane == 0) phaseObservations[row * 12u + 3] = 1;\n}\n#endif')
 out=D/(stem+'.slang');out.write_text(text);record['outputs'][str(out)]={**seal(out),'row_words':words,'clock_slots':slots}
for p in (R/'src/sirius/kernels').glob('retained_*.slang'):
 record['input_seals'][str(p)]=seal(p)
record['input_seals'][str(Path(__file__))]=seal(Path(__file__))
(W/'phase-preparation.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record['outputs']))
