"""Connect the admitted layer module to the reviewed diagnostic test observer."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, re, struct, subprocess
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
HEAD = '8f1118dd56d10cdd0f932613427b68409a1f77f6'
def seal(p):
    return {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
def one(text, old, new):
    assert text.count(old) == 1, (old, text.count(old))
    return text.replace(old, new)
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == HEAD
assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
assert json.loads((WORK / 'compiled-review.json').read_text())['accepted']
source = ROOT / 'tests/backend/shader_clock_observer.h'
before = source.read_text()
assert seal(source) == {'bytes': 643479, 'sha256': '410e7b61a77aa1007a5e0112da641725b917b2bf047cd0c426457c4e9efc0e15'}
module = WORK / 'phase/retained_transport.spv'
data = module.read_bytes()
assert seal(module) == {'bytes': 87928, 'sha256': 'e3a1589d852a4254312bb0e9e4b0343ab6aea32855c02c4f0f93b9fe1e18e5a3'}
words = struct.unpack('<' + 'I' * (len(data) // 4), data)
array = 'inline constexpr std::uint32_t kTransportClock[] = {\n'
for i in range(0, len(words), 7):
    array += '    ' + ', '.join(f'0x{v:08x}u' for v in words[i:i + 7]) + ',\n'
array += '};'
text, count = re.subn(r'inline constexpr std::uint32_t kTransportClock\[\] = \{.*?\n\};', array, before, count=1, flags=re.S)
assert count == 1
text = one(text, '// Exact admitted module SHA-256: 0eff0f81bf94582c5401e544b28bca9056ee77beb0264d4ece541aa8443dd076', '// Exact admitted module SHA-256: e3a1589d852a4254312bb0e9e4b0343ab6aea32855c02c4f0f93b9fe1e18e5a3')
text = one(text, '#include <atomic>\n', '#include <atomic>\n#include <cstring>\n')
text = one(text, 'class Observer {', '''// 16 phase cuts, followed by seven stages of 117 bounded layer records.
// Each layer retains its node range, a three-cut prefix mask and uint2 clocks.
inline constexpr std::uint32_t kTransportObservationWords = 8226;

class Observer {''')
text = one(text, 'std::uint64_t{capacity} * 36 * sizeof(std::uint32_t)', 'std::uint64_t{capacity} * kTransportObservationWords * sizeof(std::uint32_t)')
text = one(text, 'kernels_.push_back({handle, 36})', 'kernels_.push_back({handle, kTransportObservationWords})')
text = one(text, '        const auto ordinal = sequence_.fetch_add(1);', '''        // Capture the real fixed-capacity table after completion, before any
        // next command can reuse it. Active groups do not locate its program.
        // The underlying device bypasses the original probe read journal.
        std::vector<std::byte> input_snapshot;
        std::string input_sha256;
        std::uint32_t stored_capacity = 0;
        if (kernel->words == kTransportObservationWords) {
            const auto input_words = 1 + std::size_t{capacity_} * 230 +
                                     retained_program::kTransportProgram.size();
            input_snapshot.resize(input_words * sizeof(std::uint32_t));
            auto captured = device_.ReadBuffer(command.buffers[0], input_snapshot);
            if (!captured) return captured;
            std::memcpy(&stored_capacity, input_snapshot.data(), sizeof(stored_capacity));
            const auto input_digest = base::Sha256Hex(std::span(
                reinterpret_cast<const std::uint8_t*>(input_snapshot.data()),
                input_snapshot.size()));
            if (!input_digest) return Fail("cannot hash original input snapshot");
            input_sha256 = *input_digest;
        }
        const auto ordinal = sequence_.fetch_add(1);''')
text = one(text, '        const auto digest = base::Sha256Hex(', '''        if (!input_snapshot.empty()) {
            std::ofstream snapshot(directory / (stem + ".input.bin"), std::ios::binary);
            snapshot.write(reinterpret_cast<const char*>(input_snapshot.data()),
                           static_cast<std::streamsize>(input_snapshot.size()));
            snapshot.close();
            if (!snapshot) return Fail("cannot preserve complete original input snapshot");
        }
        const auto digest = base::Sha256Hex(''')
text = one(text, '''            "\\\"bytes\\\":{},\\\"sidecar_allocation_bytes\\\":{},\\\"sha256\\\":\\\"{}\\\"}}\\n",''', '''            "\\\"bytes\\\":{},\\\"sidecar_allocation_bytes\\\":{},\\\"sha256\\\":\\\"{}\\\","
            "\\\"input_bytes\\\":{},\\\"stored_capacity\\\":{},\\\"input_sha256\\\":\\\"{}\\\"}}\\n",''')
text = one(text, 'kernel->words == 36 ? "transport" : "endpoint"', 'kernel->words == kTransportObservationWords ? "transport" : "endpoint"')
text = one(text, 'sidecars_[slot].value, raw.size(), charge_, *digest);', 'sidecars_[slot].value, raw.size(), charge_, *digest, input_snapshot.size(),\n            stored_capacity, input_sha256);')
previous = WORK / 'host-before.h'
previous.write_text(before)
source.write_text(text)
record = {'parent_revision': HEAD, 'baseline_revision': '2053ef79c90752d37b0075b8137cfb9b571b6820',
          'input_seals': {str(p): seal(p) for p in (Path(__file__), module, WORK / 'compiled-review.json', previous)},
          'output': {str(source): seal(source)}, 'transport_words_per_row': 8226,
          'scope': 'Unadopted diagnostic host source only. Real sidecar charge remains inside unchanged8MiB; original tests/guards/journals/modes retained. Full actual Transport input captured using original allocated capacity; stored header must drive future table interpretation. Original Endpoint observer/module retained; host source review/build/runtime pending.'}
(WORK / 'host-preparation.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record['output']))
