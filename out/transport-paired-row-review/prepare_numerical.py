"""Prepare an isolated stateful raw/reference gate; do not execute it."""
from pathlib import Path
import hashlib
import json
import re
import shutil
import struct
import subprocess

ROOT = Path(__file__).resolve().parents[2]
TASK = Path(__file__).resolve().parent
OLD = ROOT / 'attestations/native-vulkan/3045067/transport-lanes32-rejected'
CAPACITY = 24
ROW_WORDS = 4699


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def win(path):
    return r'\\wsl.localhost\Ubuntu' + str(path.resolve()).replace('/', '\\')


revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
assert revision == '6845c8150ebf0df6992fe2609a4e9323c7f0612b'
assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT)
query = json.loads((TASK / 'native-query/query.json').read_text())
assert query['completed'] and query['dispatches'] == 0
assert json.loads((TASK / 'native-query/bridge-query.json').read_text())['all_selected_inputs_unchanged']
assert json.loads((TASK / 'native-query/terminal-query-audit.json').read_text(encoding='utf-8-sig'))['owned_births_absent']
gate = TASK / 'numerical'
gate.mkdir(exist_ok=False)
for name in ('inputs', 'observations', 'native-temp'):
    (gate / name).mkdir()
programs = json.loads((TASK / 'programs.json').read_text())
general, schwarzschild = programs['general'], programs['schwarzschild']
program = general + schwarzschild
assert (len(program), general[1], schwarzschild[1]) == (15967, 459, 459)
fixtures = {}
for family, name in [('general', 'reference_cases.json'), ('schwarzschild', 'schwarzschild_reference_cases.json')]:
    path = ROOT / 'tests/support/retained_transport' / name
    shutil.copy2(path, gate / name)
    fixtures[family] = json.loads(path.read_text())['cases']
assert [len(fixtures[f]) for f in ('general', 'schwarzschild')] == [15, 23]
header = (ROOT / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h').read_text()
words = [int(x) for x in re.findall(r'(\d+)u', re.search(r'kTransportShader\{\{(.*?)\}\};', header, re.S)[1])]
(gate / 'baseline.spv').write_bytes(struct.pack('<' + str(len(words)) + 'I', *words))
assert sha(gate / 'baseline.spv') == 'b52e6373f6b30aac6c7fe877f6e07f3c47c7ac9eafa1671840eacd19960d7cb1'
shutil.copy2(TASK / 'paired-transport.spv', gate / 'candidate.spv')
assert sha(gate / 'candidate.spv') == '721be7625f8eebbf722adf3e79a89050bffe9bf9b7260b8f5120b0a85ccd7f27'


def fixture(family, index):
    return {'words': fixtures[family][index]['input'][:], 'family': family, 'index': index,
            'name': fixtures[family][index]['name'], 'expected': 'valid'}


def changed(row, value, component, word, name, expected='refused'):
    row = {**row, 'words': row['words'][:]}
    row['words'][value * 5 + component] = word
    row.update(name=name, expected=expected)
    return row


g, s, flat = fixture('general', 0), fixture('schwarzschild', 0), fixture('general', 14)
bad = changed(s, 43, 4, 0, 'invalid-last-column')
bad_chart = changed(g, 44, 0, 0, 'invalid-chart')
bad_radius = changed(s, 0, 3, 0x80000000, 'negative-zero-radius')
# Keep first-stage products within the original2^60 factor ceiling.
# With p_x=16, next-stage x lies in (3,4)*2^60 and the RHS scale is2^62.
late = changed(s, 45, 0, 0x5d800000, 'late-arithmetic-refusal', 'late-refused')
late['words'][9 * 5] = 0x41800000
late['words'][45 * 5 + 1:45 * 5 + 4] = [0, 0, 0]
normalized = fixture('schwarzschild', 0)
normalized['words'][5:10] = [0x3f000000, 0xbf000000, 0x80000000, 0, 1]
normalized['name'] = 'normalized-spin-zero'
mixed = [g, s, flat, bad, fixture('schwarzschild', 1), fixture('general', 10),
         fixture('schwarzschild', 22), bad_chart, fixture('schwarzschild', 2), late,
         normalized, bad_radius]
mixed += [fixture('general', i) for i in range(11, 14)]
mixed += [fixture('schwarzschild', i) for i in range(16, 22)]
mixed += [bad, flat, s]
assert len(mixed) == CAPACITY
initial = [CAPACITY] + [word for row in mixed for word in row['words']] + program
(gate / 'initial.bin').write_bytes(struct.pack('<' + str(len(initial)) + 'I', *initial))
sequences = []


def sequence(name, rows, table=None):
    assert 0 < len(rows) <= CAPACITY
    upload = [CAPACITY] + [word for row in rows for word in row['words']]
    if table is not None:
        upload += [0] * ((CAPACITY - len(rows)) * 230)
        upload += table
    path = gate / 'inputs' / (name + '.bin')
    path.write_bytes(struct.pack('<' + str(len(upload)) + 'I', *upload))
    sequences.append({'name': name, 'input': win(path), 'input_sha256': sha(path),
                      'active_rows': len(rows), 'upload_bytes': path.stat().st_size,
                      'full_table_upload': table is not None,
                      'rows': [{k: v for k, v in row.items() if k != 'words'} for row in rows]})


sequence('general-fifteen', [fixture('general', i) for i in range(15)])
sequence('schwarzschild-twenty-three', [fixture('schwarzschild', i) for i in range(23)])
sequence('mixed-full', mixed)
sequence('reordered-six', [mixed[i] for i in (11, 2, 1, 5, 8, 10)])
sequence('singleton-general', [g])
sequence('singleton-schwarzschild', [s])
sequence('singleton-flat', [flat])
sequence('singleton-refused', [bad])
sequence('reordered-five', [mixed[i] for i in (5, 1, 7, 2, 8)])
sequence('mixed-three', [g, bad, s])
sequence('mixed-full-restored', mixed)
partners = [g, s, flat, bad, s, g, normalized, late]
layer_base = 43 + 5 * general[0]
assert general[layer_base:layer_base + 3] == [0, 40, 47]
table = program[:]
table[layer_base + 1] = 46
sequence('valid-mixed-layer', partners, table)


def renamed(plan, old=0, new=458):
    # Injective alpha-renaming to an unused slot preserves the ordered graph.
    used = set(plan[3:43])
    for node in range(plan[0]):
        at = 43 + node * 5
        op, destination, a, b, c = plan[at:at + 5]
        used.add(destination)
        if op >= 2:
            used.add(a)
        if op in (2, 3, 4, 5, 9, 10, 11):
            used.add(b)
        if op in (10, 11):
            used.add(c)
    assert old in used and new not in used and new < plan[1]
    result = plan[:]
    for i in range(3, 43):
        if result[i] == old:
            result[i] = new
    for node in range(plan[0]):
        at = 43 + node * 5
        op = plan[at]
        positions = [at + 1]
        if op >= 2:
            positions.append(at + 2)
        if op in (2, 3, 4, 5, 9, 10, 11):
            positions.append(at + 3)
        if op in (10, 11):
            positions.append(at + 4)
        for i in positions:
            if result[i] == old:
                result[i] = new
    return result


sequence('valid-both-slot458', partners, renamed(general) + renamed(schwarzschild))
sequence('canonical-recovery', partners, program)


def refused_partners(family):
    return [{**row, 'expected': 'table-refused'} if row['family'] == family and
            row['expected'] in ('valid', 'late-refused') and row['name'] not in ('analytic-flat', 'frozen-analytic-flat')
            else row for row in partners]


table = program[:]
assert table[43 + 5 * 40] == 0
table[43 + 5 * 40] = 12
sequence('general-opcode-refusal', refused_partners('general'), table)
sequence('general-opcode-recovery', partners, program)
table = program[:]
table[len(general) + 1] ^= 1
sequence('schwarzschild-header-refusal', refused_partners('schwarzschild'), table)
sequence('schwarzschild-header-recovery', partners, program)
table = program[:]
table[3] = general[1]
sequence('general-root-refusal', refused_partners('general'), table)
sequence('general-root-recovery', partners, program)
table = program[:]
table[len(general) + 3] = schwarzschild[1]
sequence('schwarzschild-root-refusal', refused_partners('schwarzschild'), table)
sequence('schwarzschild-root-recovery', partners, program)
table = program[:]
table[layer_base + 1] = 65
sequence('general-layer-refusal', refused_partners('general'), table)
sequence('general-layer-recovery', partners, program)
table = program[:]
table[len(general) + 43 + 5 * schwarzschild[0] + 1] = 65
sequence('schwarzschild-layer-refusal', refused_partners('schwarzschild'), table)
sequence('schwarzschild-layer-recovery', partners, program)
sequence('final-singleton-recovery', [s])
(gate / 'sequences.json').write_text(json.dumps(sequences, indent=2) + '\n')

source = (OLD / 'native_probe.cs').read_text()
assert source.count('public static class NativeTransportFmaProbe') == 1
source = source.replace('public static class NativeTransportFmaProbe', 'public static class NativePairedTransportGate')
source = source.replace('Finite native Transport comparison using the existing private Win64 Vulkan harness.',
                        'Finite stateful raw/reference gate for isolated native paired Transport.')
source = source.replace('apiVersion = (1u << 22) | (2u << 12)', 'apiVersion = (1u << 22) | (3u << 12)')
source = source.replace('canary[i] = (byte)(expected[i] ^ 255);', 'canary[i] = 255;')
needle = 'string moduleBinding, string inputBinding, bool wide, string alternativePath, string alternativeBinding) {'
assert source.count(needle) == 1
source = source.replace(needle, '''string moduleBinding, string inputBinding, bool wide, string alternativePath, string alternativeBinding,
                             string[] sequencePaths, string[] sequenceBindings, uint[] activeCounts,
                             int[] uploadBytes, string[] sequenceNames) {''')
source = source.replace('if (wide && U32(supported, 156) != 1)', 'if (U32(supported, 156) != 1)')
source = source.replace('if (wide) Marshal.WriteInt32(enabled, 156, 1);', 'Marshal.WriteInt32(enabled, 156, 1);')
source = source.replace('// Enable only the conservative binary64 product-mode feature when requested.',
                        '// Match production logical-device Float64/FMA32; modules remain native-default.')
start = source.index('            uint validBits = U32(familyList,')
end = source.index('            IntPtr priority = Allocate', start)
source = source[:start] + source[end:]
start = source.index('            var createQueryPool = Bind<CreateObject>')
end = source.index('            var resetCommand = Bind<ResetCommand>', start)
source = source[:start] + source[end:]
start = source.index('            ulong queryPool;')
end = source.index('            IntPtr queue;', start)
source = source[:start] + source[end:]
start = source.index('            var samples = new List<object>();')
end = source.index('            var loadedModules = new List<object>();', start)
source = source[:start] + '''            var samples = new List<object>(); byte[] actual = new byte[expected.Length];
            byte[] inputShadow = (byte[])input.Clone();
            if (sequencePaths.Length == 0 || sequencePaths.Length > 32 ||
                sequenceBindings.Length != sequencePaths.Length || activeCounts.Length != sequencePaths.Length ||
                uploadBytes.Length != sequencePaths.Length || sequenceNames.Length != sequencePaths.Length ||
                outputWords != checked((int)Rows * 4699)) throw new Exception("Sequence shape mismatch");
            uint dispatches = 0;
            for (int sequence = 0; sequence < sequencePaths.Length; ++sequence) {
                uint active = activeCounts[sequence];
                if (active == 0 || active > Rows ||
                    !System.Text.RegularExpressions.Regex.IsMatch(sequenceNames[sequence], "^[a-z0-9-]+$"))
                    throw new Exception("Active prefix/name invalid");
                byte[] update = ReadBounded(sequencePaths[sequence]);
                if (Hash(update) != sequenceBindings[sequence] || update.Length != uploadBytes[sequence] ||
                    Word(update, 0) != Rows ||
                    (update.Length != checked(4 + (int)active * 920) && update.Length != input.Length))
                    throw new Exception("Exact uploaded prefix/table binding mismatch");
                Array.Copy(update, 0, inputShadow, 0, update.Length);
                IntPtr uploaded; Check(mapMemory(device, memories[0], 0, (ulong)update.Length, 0, out uploaded), "vkMapMemory prefix");
                try { Marshal.Copy(update, 0, uploaded, update.Length); }
                finally { unmapMemory(device, memories[0]); }
                byte[] baseline = null;
                for (int variant = 0; variant < 2; ++variant) {
                    bool candidate = variant == 1;
                    IntPtr clear; Check(mapMemory(device, memories[1], 0, (ulong)actual.Length, 0, out clear), "vkMapMemory canary");
                    try { Marshal.Copy(canary, 0, clear, canary.Length); }
                    finally { unmapMemory(device, memories[1]); }
                    Check(resetCommand(command, 0), "vkResetCommandBuffer");
                    Check(resetFences(device, 1, fences), "vkResetFences");
                    var beginInfo = new CommandBeginInfo { sType = 42, flags = 1 };
                    Check(begin(command, ref beginInfo), "vkBeginCommandBuffer");
                    var before = new MemoryBarrier { sType = 46, source = 0x40 | 0x4000, destination = 0x20 | 0x40 };
                    barrier(command, 0x800 | 0x4000, 0x800, 0, 1, ref before, 0, IntPtr.Zero, 0, IntPtr.Zero);
                    bindPipeline(command, 1, candidate ? alternativePipeline : pipeline);
                    bindSets(command, 1, pipelineLayout, 0, 1, set, 0, IntPtr.Zero);
                    uint groupsX = candidate ? (active + 1) / 2 : active;
                    uint groupsY = candidate && (active & 1) != 0 ? 2u : 1u;
                    dispatch(command, groupsX, groupsY, 1);
                    var after = new MemoryBarrier { sType = 46, source = 0x40, destination = 0x20 | 0x2000 };
                    barrier(command, 0x800, 0x800 | 0x4000, 0, 1, ref after, 0, IntPtr.Zero, 0, IntPtr.Zero);
                    Check(end(command), "vkEndCommandBuffer");
                    submitted = false; completed = false;
                    Check(submit(queue, 1, ref submitInfo, fence), "vkQueueSubmit"); submitted = true;
                    Check(wait(device, 1, fences, 1, FenceTimeoutNanoseconds), "bounded vkWaitForFences");
                    completed = true; ++dispatches;
                    IntPtr read; Check(mapMemory(device, memories[1], 0, (ulong)actual.Length, 0, out read), "vkMapMemory read");
                    try { Marshal.Copy(read, actual, 0, actual.Length); }
                    finally { unmapMemory(device, memories[1]); }
                    byte[] inputAfter = new byte[inputShadow.Length];
                    Check(mapMemory(device, memories[0], 0, (ulong)inputAfter.Length, 0, out read), "vkMapMemory input seal");
                    try { Marshal.Copy(read, inputAfter, 0, inputAfter.Length); }
                    finally { unmapMemory(device, memories[0]); }
                    if (Hash(inputAfter) != Hash(inputShadow)) throw new Exception("Device changed input/table words");
                    string path = actualPath + "." + sequenceNames[sequence] + (candidate ? ".candidate.bin" : ".baseline.bin");
                    File.WriteAllBytes(path, actual);
                    if (baseline == null) baseline = (byte[])actual.Clone();
                    uint differences = 0, inactiveChanges = 0;
                    var examples = new List<object>();
                    for (int i = 0; i < outputWords; ++i) {
                        uint a = Word(actual, i * 4), e = Word(baseline, i * 4);
                        if (a != e) { ++differences; if (examples.Count < 8) examples.Add(Record("word", i, "expected", e, "actual", a)); }
                        if ((uint)i >= active * 4699 && a != 0xffffffffu) ++inactiveChanges;
                    }
                    samples.Add(Record("sequence", sequenceNames[sequence], "variant", candidate ? "candidate" : "baseline",
                        "active_rows", active, "groups_x", groupsX, "groups_y", groupsY,
                        "uploaded_bytes", update.Length, "input_shadow_sha256", Hash(inputShadow),
                        "output_file", Path.GetFileName(path), "output_sha256", Hash(actual),
                        "whole_capacity_word_mismatches", differences, "inactive_word_changes", inactiveChanges));
                    checkpoint("sequence=" + sequenceNames[sequence] + ";variant=" + variant + ";differences=" + differences + ";inactive=" + inactiveChanges);
                    File.AppendAllText(actualPath + ".observations.tsv", String.Join("\\t", new string[] {
                        sequenceNames[sequence], variant.ToString(), active.ToString(), groupsX.ToString(), groupsY.ToString(),
                        update.Length.ToString(), Hash(inputShadow), Hash(actual), differences.ToString(), inactiveChanges.ToString()
                    }) + Environment.NewLine);
                    if (differences != 0 || inactiveChanges != 0) throw new Exception("Paired prefix changed complete output or inactive canary; examples=" + examples.Count);
                }
            }
''' + source[end:]
start = source.index('            return Record("scope",')
end = source.index('        } finally {', start)
source = source[:start] + '''            return Record("completed", true, "scope", "Isolated paired-row raw/reference gate; no timing or full science/release acceptance",
                "device", deviceName, "driver", driverInfo, "shader_fma32_enabled", true, "shader_float64_enabled", true,
                "loaded_modules", loadedModules, "allocations", allocations, "actual_resident_bytes", allocationBytes,
                "module_sha256", moduleHash, "alternative_module_sha256", Hash(alternative),
                "capacity", Rows, "completed_dispatches", dispatches, "fence_timeout_ns", FenceTimeoutNanoseconds,
                "observed_words_per_dispatch", outputWords, "samples", samples);
''' + source[end:]
(gate / 'native_gate.cs').write_text(source)
bindings = [{'path': str(p.resolve()), 'bytes': p.stat().st_size, 'sha256': sha(p)}
            for p in [Path(__file__).resolve(), OLD / 'native_probe.cs',
                      ROOT / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h',
                      ROOT / 'tests/support/retained_transport/reference_cases.json',
                      ROOT / 'tests/support/retained_transport/schwarzschild_reference_cases.json']]
record = {'revision': revision, 'capacity': CAPACITY, 'output_words': CAPACITY * ROW_WORDS,
          'sequences': len(sequences), 'expected_dispatches': len(sequences) * 2,
          'initial_input': win(gate / 'initial.bin'), 'initial_input_sha256': sha(gate / 'initial.bin'),
          'baseline_module_sha256': sha(gate / 'baseline.spv'), 'candidate_module_sha256': sha(gate / 'candidate.spv'),
          'source_bindings': bindings, 'production_mutation': False,
          'device_dispatches_executed': 0,
          'scope': 'Prepared only. Full-capacity raw identity and inactive-canary checks on actual shared buffers; independent unchanged numerical references still required. No throughput or qualification verdict.'}
(gate / 'preparation.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'prepared': str(gate), 'sequences': len(sequences), 'expected_dispatches': len(sequences) * 2,
                  'source_bytes': len(source.encode()), 'native_execution': False}, indent=2))
