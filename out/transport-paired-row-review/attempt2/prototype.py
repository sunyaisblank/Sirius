"""Build an isolated two-row native Transport source/module; never edit production."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
TASK = Path(__file__).resolve().parent

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def identity(path):
    data = path.read_bytes()
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}

SCRATCH = '''[[vk::binding(1, 0)]] RWStructuredBuffer<uint> outputs;
// Two independent row partitions share one workgroup. Five words remain
// explicit. The complete 459-register domain includes an output-backed spill.
groupshared uint retainedScratch[2 * 384 * 5];
groupshared uint retainedStatus[2];
uint RetainedPartition() {
    return spirv_asm { OpLoad $$uint result builtin(LocalInvocationIndex:uint); } / 32;
}
uint ScratchRead(uint index, uint component, uint rowBase, uint prefix) {
    if (index < 384)
        return retainedScratch[(RetainedPartition() * 384 + index) * 5 + component];
    return outputs[rowBase + prefix + index * 5 + component];
}
void ScratchWrite(uint index, uint component, uint rowBase, uint prefix, uint value) {
    if (index < 384)
        retainedScratch[(RetainedPartition() * 384 + index) * 5 + component] = value;
    else
        outputs[rowBase + prefix + index * 5 + component] = value;
}
void RejectRetainedProgram() {
    uint previous;
    InterlockedOr(retainedStatus[RetainedPartition()], 1, previous);
}
bool RetainedFailed() { return (retainedStatus[RetainedPartition()] & 1) != 0; }
bool AnyRetainedLive() {
    return (retainedStatus[0] & 0x80000001u) == 0x80000000u ||
           (retainedStatus[1] & 0x80000001u) == 0x80000000u;
}
'''

EVALUATE = '''bool EvaluateProgram(uint row, uint rowBase, uint programBase, uint count, uint registers,
                     uint layers, uint stage, uint lane, bool alive) {
    uint layerBase = programBase + 43 + 5 * count;
    if (alive && (count == 0 || layers == 0 || layers > SIRIUS_RETAINED_LAYERS ||
                  inputs[layerBase] != 0 || inputs[layerBase + layers] != count))
        RejectRetainedProgram();
    AllMemoryBarrierWithGroupSync();
    alive = alive && !RetainedFailed();
    uint expected = 0;
    bool done = false;
    // Both row partitions reach every collective, even after one has refused
    // or its shorter table has ended. Only live row instructions execute.
    for (uint layer = 0; layer < SIRIUS_RETAINED_LAYERS; ++layer) {
        uint start = 0, end = 0;
        if (alive && !done && layer < layers) {
            start = inputs[layerBase + layer];
            end = inputs[layerBase + layer + 1];
            if (start >= count) done = true;
            else if (start != expected || end <= start || end > count || end - start > 64)
                RejectRetainedProgram();
            else
                for (uint node = start + lane; node < end; node += 32)
                    if (!EvaluateOne(node, row, rowBase, programBase, registers))
                        RejectRetainedProgram();
        }
        // The spill is StorageBuffer memory; a workgroup-only fence is insufficient.
        AllMemoryBarrierWithGroupSync();
        alive = alive && !RetainedFailed();
        if (alive && !done && layer < layers) expected = end;
    }
    if (alive)
        for (uint i = lane; i < 40; i += 32) {
            uint index = inputs[programBase + 3 + i];
            if (index >= registers) RejectRetainedProgram();
            else if (ScratchRead(index, 4, rowBase, 2404) == 0) RejectRetainedProgram();
        }
    AllMemoryBarrierWithGroupSync();
    alive = alive && !RetainedFailed();
    if (alive)
        for (uint i = lane; i < 40; i += 32)
            WriteStored(rowBase + 1004 + stage * 200 + i * 5,
                        ReadValue(inputs[programBase + 3 + i], rowBase));
    AllMemoryBarrierWithGroupSync();
    return alive;
}
'''

MAIN = '''bool WriteFlatRow(uint row, uint rowBase, RetainedTriple h) {
    for (uint i = 0; i < 40; ++i) {
        RetainedTriple delta = RTConstant(RPZero());
        if (i % 8 < 4) {
            RetainedTriple slope = ReadOriginal(8 + i, row);
            if (i % 8 == 0) slope = RTNegate(slope);
            delta = RTMultiply(h, slope);
        }
        RetainedTriple full = RTAdd(ReadOriginal(4 + i, row), delta);
        if (full.valid == 0 || delta.valid == 0) return false;
        WriteStored(rowBase + 4 + i * 5, full);
        WriteStored(rowBase + 204 + i * 5, full);
        WriteStored(rowBase + 404 + i * 5, delta);
        WriteStored(rowBase + 604 + i * 5, RTConstant(RPZero()));
    }
    outputs[rowBase] = 1;
    outputs[rowBase + 1] = 7;
    return true;
}
[shader("compute")][numthreads(64, 1, 1)] void ComputeMain(uint3 id : SV_GroupID,
                                                        uint globalLane : SV_GroupIndex) {
    uint3 dispatch = spirv_asm { OpLoad $$uint3 result builtin(NumWorkgroups:uint3); };
    uint ni, si, no, so;
    inputs.GetDimensions(ni, si);
    outputs.GetDimensions(no, so);
    // y=2 encodes an odd active prefix. Its y=1 workgroups never touch memory.
    if (ni < 1 || id.y != 0 || id.z != 0 || dispatch.z != 1 ||
        (dispatch.y != 1 && dispatch.y != 2)) return;
    uint rows = inputs[0], row = id.x * 2 + globalLane / 32, lane = globalLane % 32;
    uint activeRows = min(rows, dispatch.x * 2 - (dispatch.y == 2 ? 1 : 0));
    if (rows == 0 || rows > 65536 || id.x * 2 >= activeRows ||
        ni < 1 + 230 * rows + 43 || no % rows != 0) return;
    bool real = row < activeRows;
    uint programBase = 1 + 230 * rows, count = inputs[programBase],
         registers = inputs[programBase + 1];
    uint rowWords = no / rows, rowBase = row * rowWords;
    if (rowWords < 2404) return;
    if (real)
        for (uint i = lane; i < rowWords; i += 32) outputs[rowBase + i] = 0;
    AllMemoryBarrierWithGroupSync();
    uint layers = SIRIUS_RETAINED_LAYERS;
    if (inputs[programBase + 2] != 40 || count > 65536 || registers != SIRIUS_RETAINED_REGISTERS ||
        43 + 5 * count + layers + 1 != SIRIUS_RETAINED_GENERAL_PROGRAM_WORDS ||
        ni != programBase + SIRIUS_RETAINED_GENERAL_PROGRAM_WORDS +
                            SIRIUS_RETAINED_SCHWARZSCHILD_PROGRAM_WORDS ||
        rowWords != 2404 + 5 * registers) return;
    bool alive = real;
    for (uint i = 0; alive && i < 46; ++i)
        if (ReadOriginal(i, row).valid == 0) alive = false;
    RetainedTriple h = RTInvalid(), spin = RTInvalid(), lambda = RTInvalid();
    if (alive) {
        RetainedTriple chart = ReadOriginal(44, row);
        alive = (RPScalarEqual(chart.hi, RPOne()) || RPScalarEqual(chart.hi, RPScalarNegate(RPOne()))) &&
                RPScalarEqual(chart.lo, RPZero()) && RPScalarEqual(chart.tail, RPZero()) &&
                RPScalarEqual(chart.error, RPZero());
    }
    if (alive) {
        spin = ReadOriginal(1, row); lambda = ReadOriginal(3, row);
        alive = RPScalarEqual(RTAbsoluteUpper(spin), RPZero()) ||
                RPScalarEqual(RTAbsoluteUpper(lambda), RPZero());
    }
    if (alive) {
        h = ReadOriginal(45, row);
        alive = RPScalarGreater(h.hi, RPZero()) && RPScalarGreater(RTAbsoluteLower(h), RPZero());
    }
    bool flat = false;
    if (alive) flat = RTExact(ReadOriginal(0, row), RPZero()) &&
                      RTExact(ReadOriginal(2, row), RPZero()) && RTExact(lambda, RPZero());
    if (flat && lane == 0) WriteFlatRow(row, rowBase, h);
    alive = alive && !flat;
    bool schwarzschild = false;
    if (alive) schwarzschild = RTExact(spin, RPZero()) && RTExact(ReadOriginal(2, row), RPZero()) &&
                               RTExact(lambda, RPZero());
    if (schwarzschild) {
        programBase += SIRIUS_RETAINED_GENERAL_PROGRAM_WORDS;
        count = inputs[programBase]; layers = SIRIUS_RETAINED_SCHWARZSCHILD_LAYERS;
        if (inputs[programBase + 1] != registers || inputs[programBase + 2] != 40 ||
            count > 65536 || 43 + 5 * count + layers + 1 != SIRIUS_RETAINED_SCHWARZSCHILD_PROGRAM_WORDS)
            alive = false;
    }
    if (lane == 0) retainedStatus[globalLane / 32] = alive ? 0x80000000u : 0;
    if (alive)
        for (uint index = lane; index < SIRIUS_RETAINED_REGISTERS; index += 32)
            ScratchWrite(index, 4, rowBase, 2404, 0);
    AllMemoryBarrierWithGroupSync();
    for (uint stage = 0; stage < 7; ++stage) {
        if (!AnyRetainedLive()) break; // uniform across the complete workgroup
        if (lane == 0) retainedStatus[globalLane / 32] = alive ? 0x80000000u : 0;
        AllMemoryBarrierWithGroupSync();
        // One workgroup producer per coefficient, preserving the original value.
        if (stage == 1) CacheCoefficient(globalLane);
        AllMemoryBarrierWithGroupSync();
        if (alive)
            for (uint i = lane; i < 40; i += 32) {
                RetainedTriple value = ReadOriginal(4 + i, row);
                if (stage != 0) {
                    RetainedTriple first = ReadStored(rowBase + 1004 + i * 5);
                    RetainedTriple weighted = RTMultiply(first, ReadCoefficient(stage - 1));
                    for (uint j = 1; j < stage; ++j)
                        if (an[stage][j] != 0) {
                            RetainedTriple difference =
                                RTSubtract(ReadStored(rowBase + 1004 + j * 200 + i * 5), first);
                            uint index = 6 + (stage - 2) * (stage - 1) / 2 + j - 1;
                            if (stage == 6) --index;
                            weighted = RTAdd(weighted, RTMultiply(difference, ReadCoefficient(index)));
                        }
                    RetainedTriple increment = RTMultiply(weighted, h);
                    value = RTAdd(value, increment);
                    if (stage == 6) WriteStored(rowBase + 4 + 80 * 5 + i * 5, increment);
                }
                if (value.valid == 0) RejectRetainedProgram();
                WriteStored(rowBase + 804 + i * 5, value);
            }
        AllMemoryBarrierWithGroupSync();
        alive = alive && !RetainedFailed();
        if (alive && lane == 0) outputs[rowBase + 1] = stage + 1;
        alive = EvaluateProgram(row, rowBase, programBase, count, registers, layers, stage, lane, alive);
        if (lane == 0) retainedStatus[globalLane / 32] = alive ? 0x80000000u : 0;
        AllMemoryBarrierWithGroupSync();
    }
    if (lane == 0) retainedStatus[globalLane / 32] = alive ? 0x80000000u : 0;
    AllMemoryBarrierWithGroupSync();
    if (alive)
        for (uint i = lane; i < 40; i += 32) {
            RetainedTriple first = ReadStored(rowBase + 1004 + i * 5), weighted = RTConstant(RPZero());
            for (uint j = 1; j < 7; ++j)
                if (en[j] != 0) {
                    RetainedTriple difference =
                        RTSubtract(ReadStored(rowBase + 1004 + j * 200 + i * 5), first);
                    weighted = RTAdd(weighted, RTMultiply(difference, ReadCoefficient(18 + j)));
                }
            RetainedTriple error = RTMultiply(weighted, h), full = ReadStored(rowBase + 804 + i * 5);
            RetainedTriple lower = RTSubtract(full, error);
            if (error.valid == 0 || lower.valid == 0) RejectRetainedProgram();
            WriteStored(rowBase + 4 + i * 5, full);
            WriteStored(rowBase + 4 + 40 * 5 + i * 5, lower);
            WriteStored(rowBase + 4 + 120 * 5 + i * 5, error);
        }
    AllMemoryBarrierWithGroupSync();
    alive = alive && !RetainedFailed();
    if (alive && lane == 0) { outputs[rowBase] = 1; outputs[rowBase + 2] = 1; }
    // Restore the original reserved zero words for every real success/refusal row.
    if (real)
        for (uint i = lane; i < (SIRIUS_RETAINED_REGISTERS - 384) * 5; i += 32)
            outputs[rowBase + 2404 + 384 * 5 + i] = 0;
    AllMemoryBarrierWithGroupSync();
}
'''

def main():
    kernels = TASK / 'kernels'
    kernels.mkdir(exist_ok=False)
    bindings = []
    for path in sorted((ROOT / 'src/sirius/kernels').iterdir()):
        if path.is_file() and path.suffix in ('.slang', '.h', '.py'):
            shutil.copyfile(path, kernels / path.name)
            bindings.append(identity(path))
    transport = (kernels / 'retained_transport.slang').read_text()
    start = transport.index('bool EvaluateProgram(')
    end = transport.index('RetainedTriple Rational(', start)
    entry = transport.index('[shader("compute")][numthreads(')
    original_one = transport[transport.index('bool EvaluateOne('):start]
    tables = transport[end:entry]
    candidate = transport[:start] + EVALUATE + tables + MAIN
    declaration = '[[vk::binding(1, 0)]] RWStructuredBuffer<uint> outputs;\n'
    assert candidate.count(declaration) == 1
    candidate = candidate.replace(declaration, '')
    assert original_one in candidate and tables in candidate
    (kernels / 'retained_transport.slang').write_text(candidate)
    (kernels / 'retained_scratch.slang').write_text(SCRATCH)
    program = load('retained_program', ROOT / 'src/sirius/kernels/retained_program.py')
    builder = load('retained_builder', ROOT / 'scripts/build-retained-kernels.py')
    general = program.build_transport_program(parallel=True)
    schwarzschild = program.build_schwarzschild_transport_program(parallel=True)
    def encode(p):
        return [p['instructions'], p['registers'], 40] + p['outputs'] + p['operations'] + p['layer_offsets']
    g, s = encode(general), encode(schwarzschild)
    (TASK / 'programs.json').write_text(json.dumps({'general': g, 'schwarzschild': s}))
    cache = (ROOT / 'bin/linux-gcc/CMakeCache.txt').read_text()
    def tool(name):
        return re.search('^SIRIUS_' + name + ':FILEPATH=(.*)$', cache, re.M)[1]
    paths = {n: tool(n) for n in ('SLANGC', 'SPIRV_AS', 'SPIRV_DIS', 'SPIRV_VAL')}
    paths['SPIRV_OPT'] = re.search('^SIRIUS_SPIRV_OPT:FILEPATH=(.*)$', cache, re.M)[1]
    bindings += [identity(Path(p)) for p in paths.values()]
    bindings += [identity(ROOT / 'scripts/build-retained-kernels.py'), identity(Path(__file__))]
    receipt = {'baseline': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
               'git_status': subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True),
               'inputs': bindings, 'prototype_only': True, 'device_dispatches': 0,
               'logical_shared_bytes': 2 * 384 * 5 * 4 + 25 * 5 * 4 + 2 * 4,
               'unchanged_evaluate_one_sha256': hashlib.sha256(original_one.encode()).hexdigest(),
               'unchanged_coefficients_sha256': hashlib.sha256(tables.encode()).hexdigest(),
               'candidate_sources': [identity(kernels / name) for name in ('retained_transport.slang', 'retained_scratch.slang')]}
    (TASK / 'compile-start.json').write_text(json.dumps(receipt, indent=2) + '\n')
    started = time.monotonic()
    module = TASK / 'paired-transport.spv'
    try:
        builder.compile_shader(kernels / 'retained_transport.slang', module,
            paths['SLANGC'], paths['SPIRV_AS'], paths['SPIRV_DIS'], paths['SPIRV_VAL'],
            general['registers'], 5, len(general['layer_offsets']) - 1,
            optimizer=paths['SPIRV_OPT'], transport_specialized=(len(g), len(s), len(schwarzschild['layer_offsets']) - 1))
        subprocess.run([paths['SPIRV_DIS'], str(module), '-o', str(TASK / 'paired-transport.spvasm')], check=True)
        receipt['module'] = identity(module)
        receipt['assembly'] = identity(TASK / 'paired-transport.spvasm')
        receipt['compile_success'] = True
    except Exception as error:
        receipt['compile_success'] = False
        receipt['error'] = repr(error)
        raise
    finally:
        receipt['wall_seconds'] = time.monotonic() - started
        receipt['all_selected_inputs_unchanged'] = all(identity(Path(record['path'])) == record for record in bindings)
        (TASK / 'compile-result.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ('compile_success', 'wall_seconds', 'logical_shared_bytes', 'all_selected_inputs_unchanged', 'module')}, indent=2), flush=True)

if __name__ == '__main__':
    main()
