"""Prepare test-only word witnesses from reviewed modules; no build/device."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib
import json
import struct
import subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
BASE = '9fbd0e8b0d549a13476d066f2877dcd2cf9a94d9'
OLD_BASE = '2053ef79c90752d37b0075b8137cfb9b571b6820'
OLD_DIAGNOSTIC = '24a8e145d1da5dca75a4fa26ec8e639020eecf37'


def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def one(text, old, new):
    assert text.count(old) == 1, (old, text.count(old))
    return text.replace(old, new)


assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == BASE
assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
compiled_review = json.loads((WORK / 'compiled-review.json').read_text())
assert compiled_review['accepted']
assert compiled_review['source_revision'] == BASE
assert compiled_review['schema_seal'] == seal(WORK / 'value-schema.json')
for stage in ('transport', 'endpoint'):
    assert compiled_review['modules'][stage] == seal(WORK / 'kernels' / ('retained_' + stage + '.spv'))
old_header = ROOT / 'attestations/native-vulkan/24a8e14/native-rhs-layer-clock-execution/files/host-before-format.h'
old = old_header.read_text()
prefix = old[:old.index('// Exact admitted module SHA-256:')]
prefix = prefix.replace('// Temporary diagnostic: original native-default modules with lane-zero clock\n// writes. No arithmetic, output, refusal, barrier, or governor changes.',
                        '// Private diagnostic: original native-default modules with word snapshots.\n// Numerical operations, original outputs, guards and barriers are retained.')
prefix = prefix.replace('shader_clock', 'value_witness')
prefix = prefix.replace('#include <span>\n', '#include <span>\n#include <string>\n')
modules = {}
arrays = ''
for stage in ('transport', 'endpoint'):
    module = WORK / 'kernels' / ('retained_' + stage + '.spv')
    data = module.read_bytes()
    assert len(data) % 4 == 0
    words = struct.unpack('<' + 'I' * (len(data) // 4), data)
    modules[stage] = seal(module)
    name = 'k' + stage.title() + 'Witness'
    arrays += '// Exact observed module SHA-256: ' + modules[stage]['sha256'] + '\n'
    arrays += 'inline constexpr std::uint32_t ' + name + '[] = {\n'
    for offset in range(0, len(words), 7):
        arrays += '    ' + ', '.join(f'0x{v:08x}u' for v in words[offset:offset + 7]) + ',\n'
    arrays += '};\n'

body = old[old.index('class Observer {'):]
body = body.replace('shader_clock', 'value_witness').replace('ShaderClock', 'ValueWitness')
body = body.replace('shader-clock', 'word-value').replace('kTransportObservationWords', 'kWordsPerRow')
body = body.replace('kTransportClock', 'kTransportWitness').replace('kEndpointClock', 'kEndpointWitness')
body = one(body, 'kernels_.push_back({handle, kWordsPerRow})',
           'kernels_.push_back({handle, Stage::kTransport})')
body = one(body, 'kernels_.push_back({handle, 12})',
           'kernels_.push_back({handle, Stage::kEndpoint})')
body = one(body, '    struct Kernel {\n        KernelHandle handle;\n        std::uint32_t words;\n    };',
           '    enum class Stage { kTransport, kEndpoint };\n    struct Kernel {\n        KernelHandle handle;\n        Stage stage;\n    };')
body = one(body,
           '            "Perturbed diagnostic calls include sidecar clear/readback and input snapshot/hash/file/XML work; "\n'
           '            "underlying device markers and DispatchTiming omit host observer work. "\n'
           '            "Raw device ticks are uncalibrated; no exclusive-cost or ordinary-benefit claim.");',
           '            "Perturbed diagnostic calls include sidecar clear/readback, complete input/output snapshots, "\n'
           '            "hash/file/XML work; underlying device markers and DispatchTiming omit host observer work. "\n'
           '            "No ordinary-performance, universal-word-equivalence or qualification claim.");')
capture = r'''    base::Expected<void> Capture(const ComputeDispatch& command, std::size_t slot, bool paired) {
        const auto* kernel = Find(command.kernel);
        if (!kernel) return {};
        if (slot >= sidecars_.size() || command.buffers.size() != 2)
            return Fail("invalid completed observer binding");
        const bool transport = kernel->stage == Stage::kTransport;
        const auto input_row_words = transport ? 230U : 225U;
        const auto output_row_words = transport ? 4699U : 3829U;
        const auto program_words = transport ? retained_program::kTransportProgram.size()
                                             : retained_program::kEndpointProgram.size();
        // Synchronous dispatch has completed. Capture both complete original
        // allocations before callbacks or the next call can reuse their words.
        // These direct device reads leave the original probe journal intact.
        std::vector<std::byte> raw(zero_.size());
        std::vector<std::byte> input_snapshot(
            (1 + std::size_t{capacity_} * input_row_words + program_words) * sizeof(std::uint32_t));
        std::vector<std::byte> output_snapshot(
            std::size_t{capacity_} * output_row_words * sizeof(std::uint32_t));
        auto captured = device_.ReadBuffer(sidecars_[slot], raw);
        if (!captured) return captured;
        captured = device_.ReadBuffer(command.buffers[0], input_snapshot);
        if (!captured) return captured;
        captured = device_.ReadBuffer(command.buffers[1], output_snapshot);
        if (!captured) return captured;
        std::uint32_t stored_capacity = 0;
        std::memcpy(&stored_capacity, input_snapshot.data(), sizeof(stored_capacity));
        const auto digest = Hash(raw);
        const auto input_digest = Hash(input_snapshot);
        const auto output_digest = Hash(output_snapshot);
        if (!digest || !input_digest || !output_digest)
            return Fail("cannot hash complete diagnostic snapshots");
        const auto ordinal = sequence_.fetch_add(1);
        const auto* test = ::testing::UnitTest::GetInstance()->current_test_info();
        if (!test) return Fail("missing original test identity");
        const std::string stem = std::format("{}_{}_{}", test->test_suite_name(), test->name(), ordinal);
        const std::filesystem::path directory("word-value-observations");
        std::error_code error;
        std::filesystem::create_directories(directory, error);
        if (error) return Fail("cannot create diagnostic capture directory");
        if (std::filesystem::exists(directory / (stem + ".json"), error) || error)
            return Fail("diagnostic capture identity already exists");
        auto preserved = Preserve(directory / (stem + ".bin"), raw);
        if (!preserved) return preserved;
        preserved = Preserve(directory / (stem + ".input.bin"), input_snapshot);
        if (!preserved) return preserved;
        preserved = Preserve(directory / (stem + ".output.bin"), output_snapshot);
        if (!preserved) return preserved;
        const std::string metadata = std::format(
            "{{\"version\":1,\"ordinal\":{},\"stage\":\"{}\",\"capacity\":{},\"active_rows\":{},"
            "\"groups_y\":{},\"groups_z\":{},\"words_per_row\":{},\"slot\":{},"
            "\"paired\":{},\"kernel\":{},\"input\":{},\"output\":{},\"sidecar\":{},"
            "\"bytes\":{},\"sidecar_allocation_bytes\":{},\"sha256\":\"{}\","
            "\"input_bytes\":{},\"stored_capacity\":{},\"input_sha256\":\"{}\","
            "\"output_bytes\":{},\"output_sha256\":\"{}\",\"module_sha256\":\"{}\","
            "\"schema_sha256\":\"SCHEMA_SHA\",\"input_row_words\":{},\"output_row_words\":{}}}\n",
            ordinal, transport ? "transport" : "endpoint", capacity_, command.groups_x,
            command.groups_y, command.groups_z, kWordsPerRow, slot, paired ? "true" : "false",
            command.kernel.value, command.buffers[0].value, command.buffers[1].value,
            sidecars_[slot].value, raw.size(), charge_, *digest, input_snapshot.size(),
            stored_capacity, *input_digest, output_snapshot.size(), *output_digest,
            transport ? "TRANSPORT_SHA" : "ENDPOINT_SHA", input_row_words, output_row_words);
        std::ofstream receipt(directory / (stem + ".json"), std::ios::binary);
        receipt << metadata;
        receipt.close();
        if (!receipt) return Fail("cannot preserve diagnostic capture receipt");
        ::testing::Test::RecordProperty("value_witness_capture_" + std::to_string(ordinal), metadata);
        std::clog << "[ValueWitness] " << metadata << std::flush;
        return {};
    }

'''
capture = capture.replace('SCHEMA_SHA', seal(WORK / 'value-schema.json')['sha256'])
capture = capture.replace('TRANSPORT_SHA', modules['transport']['sha256'])
capture = capture.replace('ENDPOINT_SHA', modules['endpoint']['sha256'])
start = body.index('    base::Expected<void> Capture(')
end = body.index('  private:', start)
body = body[:start] + capture + body[end:]
helpers = r'''    static std::expected<std::string, std::string> Hash(std::span<const std::byte> bytes) {
        return base::Sha256Hex(std::span(reinterpret_cast<const std::uint8_t*>(bytes.data()), bytes.size()));
    }
    static base::Expected<void> Preserve(const std::filesystem::path& path, std::span<const std::byte> bytes) {
        std::ofstream output(path, std::ios::binary);
        output.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
        output.close();
        if (!output) return Fail("cannot preserve complete diagnostic snapshot");
        return {};
    }
'''
body = one(body, '    ComputeDevice& device_;', helpers + '    ComputeDevice& device_;')
header = prefix + arrays + '\ninline constexpr std::uint32_t kWordsPerRow = 768;\n\n' + body
assert 'Clock' not in header and 'clock' not in header
assert 'kernel->words' not in header
header_path = WORK / 'word_value_observer.h'
assert not header_path.exists()
header_path.write_text(header)

test_path = 'tests/backend/retained_compute_test.cpp'
baseline = subprocess.check_output(['git', 'show', OLD_BASE + ':' + test_path])
assert (ROOT / test_path).read_bytes() == baseline
test = subprocess.check_output(['git', 'show', OLD_DIAGNOSTIC + ':' + test_path], text=True)
test = test.replace('shader_clock_observer.h', 'word_value_observer.h')
test = test.replace('sirius::test::shader_clock::Observer clock_', 'sirius::test::value_witness::Observer words_')
test = test.replace('clock_', 'words_').replace('PhaseObserver', 'WordObserver')
test = test.replace('PhaseAllocationBytes', 'WordAllocationBytes').replace('PhaseKernelCount', 'WordKernelCount')
test = test.replace('phase_enabled', 'word_enabled')
test = test.replace('sidecar clear/read/hash/file/', 'sidecar clear/read, complete input/output snapshots/hash/file/')
test = test.replace('sidecar clear/read/hash/file/XML', 'sidecar clear/read, complete input/output snapshots/hash/file/XML')
test_output = WORK / 'retained_compute_test.cpp'
assert not test_output.exists()
test_output.write_text(test)
record = {'source_revision': BASE, 'prior_diagnostic_revision': OLD_DIAGNOSTIC,
          'input_seals': {str(p.relative_to(ROOT)): seal(p) for p in
                          (Path(__file__), old_header, WORK / 'compiled-review.json', WORK / 'value-schema.json',
                           WORK / 'kernels/retained_transport.spv', WORK / 'kernels/retained_endpoint.spv', ROOT / test_path)},
          'outputs': {p.name: seal(p) for p in (header_path, test_output)},
          'modules': modules, 'words_per_row': 768,
          'scope': 'Prepared private test-only host source. Two distinct full-capacity sidecars cleared per command; actual allocation padding charged under original8MiB. Explicit stage identity; complete real original inputs and outputs captured synchronously before callbacks. Original test inputs, modes, numerical expectations,14factoryallocations, journals and ordered markers retained. Compilation/native use pending independent review.'}
(WORK / 'host-preparation.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record['outputs'], indent=2))
