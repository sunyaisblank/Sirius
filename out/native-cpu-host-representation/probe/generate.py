"""Generate one narrow comparison from the actual old and patched lambdas."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
SOURCE = R / 'src/sirius/backend/retained_compute.cpp'
PATCH = W.parent / 'host.patch'


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


def apply(text, patch):
    lines, changes = text.splitlines(keepends=True), patch.splitlines(keepends=True)
    assert changes[:2] == ['--- a/src/sirius/backend/retained_compute.cpp\n',
                           '+++ b/src/sirius/backend/retained_compute.cpp\n']
    result, position, cursor = [], 0, 2
    while cursor < len(changes):
        match = re.fullmatch(r'@@ -(\d+),(\d+) \+(\d+),(\d+) @@.*\n', changes[cursor])
        assert match
        first, old_count, _, new_count = map(int, match.groups())
        assert first - 1 >= position
        result.extend(lines[position:first - 1]); position = first - 1
        cursor += 1
        removed = added = 0
        while cursor < len(changes) and not changes[cursor].startswith('@@'):
            line = changes[cursor]; cursor += 1
            assert line[0] in ' +-'
            if line[0] in ' -':
                assert lines[position] == line[1:]
                position += 1; removed += 1
            if line[0] in ' +':
                result.append(line[1:]); added += 1
        assert (removed, added) == (old_count, new_count)
    return ''.join(result + lines[position:])


def spacing(text):
    start = text.index('bool RetainedValue::IsRepresented() const {')
    match = re.match(r'bool RetainedValue::IsRepresented\(\) const \{\n'
                     r'(    const auto spacing = \[\]\(float value\) \{.*?\n    \};)',
                     text[start:], re.S)
    assert match
    return match[1]


CPP = r'''#include <algorithm>
#include <array>
#include <bit>
#include <cerrno>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <immintrin.h>
#include <limits>
#include <stdexcept>
#include <vector>

__attribute__((noinline, noipa)) double reference_spacing(float value) {
@REFERENCE@
    return spacing(value);
}
__attribute__((noinline, noipa)) double candidate_spacing(float value) {
@CANDIDATE@
    return spacing(value);
}

struct Observed {
    std::uint64_t bits;
    int error, flags;
    unsigned csr;
    int rounding;
    bool operator==(const Observed&) const = default;
};
struct Failure {
    std::uint32_t input;
    int rounding, seed;
    unsigned controls;
    Observed reference, candidate;
};
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
__attribute__((noinline, noipa))
Observed observe(double (*function)(float), std::uint32_t input, int rounding,
                 unsigned controls, int seed, const std::fenv_t& baseline) {
    require(std::fesetenv(&baseline) == 0 && std::fesetround(rounding) == 0,
            "cannot install environment");
    require(std::feclearexcept(FE_ALL_EXCEPT) == 0, "cannot clear flags");
    if (seed) require(std::feraiseexcept(FE_INEXACT) == 0, "cannot seed inexact");
    _mm_setcsr((_mm_getcsr() & ~(0x8040u | 0x3fu)) | controls | (seed ? 0x20u : 0u));
    require(std::fetestexcept(FE_ALL_EXCEPT) == (seed ? FE_INEXACT : 0) &&
            std::fegetround() == rounding && (_mm_getcsr() & 0x8040u) == controls,
            "environment precondition differs");
    const float value = std::bit_cast<float>(input);
    int* const error_pointer = &errno;
    *error_pointer = EDOM;  // A nonzero sentinel detects unexpected clearing.
    asm volatile("" ::: "memory");
    const double result = function(value);
    asm volatile("" ::: "memory");
    const int error = *error_pointer;
    const unsigned csr = _mm_getcsr();
    const int flags = std::fetestexcept(FE_ALL_EXCEPT);
    const int mode = std::fegetround();
    return {std::bit_cast<std::uint64_t>(result), error, flags, csr, mode};
}
void print_observed(const Observed& value) {
    std::printf("{\"bits\":\"%016llx\",\"errno\":%d,\"FEflags\":%d,\"MXCSR\":%u,\"round\":%d}",
                static_cast<unsigned long long>(value.bits), value.error, value.flags, value.csr, value.rounding);
}
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    int* const error_pointer = &errno;
    const int saved_errno = *error_pointer;
    std::fenv_t ambient, baseline;
    if (std::fegetenv(&ambient) != 0) return 2;
    const unsigned ambient_csr = _mm_getcsr();
    const int ambient_flags = std::fetestexcept(FE_ALL_EXCEPT);
    const int ambient_round = std::fegetround();
    const auto restore = [&] {
        const bool success = std::fesetenv(&ambient) == 0;
        _mm_setcsr(ambient_csr); *error_pointer = saved_errno;
        return success;
    };
    std::uint64_t comparisons = 0, mismatches = 0;
    std::size_t input_count = 0;
    unsigned mxcsr_mask = 0;
    bool trap_masked = false;
    const char* setup_error = nullptr;
    std::vector<Failure> failures;
    try {
        // Save actual mapped paths; Python later pins them against preparation.
        FILE* maps = std::fopen("/proc/self/maps", "r");
        FILE* copy = std::fopen(argv[1], "wx");
        require(maps && copy, "cannot preserve actual process maps");
        char line[4096];
        while (std::fgets(line, sizeof(line), maps)) std::fputs(line, copy);
        const bool io_ok = !std::ferror(maps) && !std::ferror(copy);
        const int close_maps = std::fclose(maps), close_copy = std::fclose(copy);
        require(io_ok && close_maps == 0 && close_copy == 0, "cannot finish process maps");
        alignas(16) unsigned char state[512]{};
        _fxsave(state); std::memcpy(&mxcsr_mask, state + 28, sizeof(mxcsr_mask));
        if (!mxcsr_mask) mxcsr_mask = 0xffbfu;
        require((mxcsr_mask & 0x8040u) == 0x8040u, "FTZ/DAZ unsupported");
        require(std::feholdexcept(&baseline) == 0, "cannot mask traps");
        trap_masked = true;
        require(std::fegetenv(&baseline) == 0, "cannot preserve masked baseline");
        std::vector<std::uint32_t> inputs;
        const auto add_signed = [&](std::uint32_t bits) {
            inputs.push_back(bits); inputs.push_back(bits | 0x80000000u);
        };
        for (unsigned exponent = 1; exponent <= 254; ++exponent)
            for (unsigned mantissa : {0u, 1u, 0x400000u, 0x7ffffeu, 0x7fffffu})
                add_signed((exponent << 23) | mantissa);
        for (unsigned bits : {0u, 1u, 2u, 0x003fffffu, 0x00400000u, 0x007ffffeu,
                              0x007fffffu, 0x00800000u, 0x00800001u, 0x00800002u,
                              0x7f7ffffeu, 0x7f7fffffu, 0x7f800000u, 0x7fc00000u,
                              0x7fc00001u, 0x7fffffffu, 0x7f800001u, 0x7f800002u,
                              0x7fa00000u, 0x7fbfffffu})
            add_signed(bits);
        std::sort(inputs.begin(), inputs.end());
        inputs.erase(std::unique(inputs.begin(), inputs.end()), inputs.end());
        input_count = inputs.size();
        for (int rounding : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO})
            for (unsigned mode = 0; mode != 4; ++mode)
                for (int seed : {0, 1}) {
                    const unsigned controls = (mode & 1 ? 0x40u : 0u) | (mode & 2 ? 0x8000u : 0u);
                    for (std::uint32_t input : inputs) {
                        const auto a = observe(reference_spacing, input, rounding, controls, seed, baseline);
                        const auto b = observe(candidate_spacing, input, rounding, controls, seed, baseline);
                        ++comparisons;
                        if (!(a == b)) {
                            ++mismatches;
                            if (failures.size() < 32) failures.push_back({input, rounding, seed, controls, a, b});
                        }
                    }
                }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "setup/probe error: %s\n", error.what());
        setup_error = "setup_or_probe_exception";
    } catch (...) {
        setup_error = "unknown_exception";
    }
    const bool restore_ok = restore();
    const bool restored = restore_ok && _mm_getcsr() == ambient_csr &&
                          std::fetestexcept(FE_ALL_EXCEPT) == ambient_flags &&
                          std::fegetround() == ambient_round && *error_pointer == saved_errno;
    const bool passed = !setup_error && comparisons == input_count * 32 && !mismatches && restored;
    std::printf("{\"passed\":%s,\"unique_inputs\":%zu,\"environments\":32,\"comparisons\":%llu,"
                "\"function_calls\":%llu,\"mismatches\":%llu,\"MXCSR_supported_mask\":%u,"
                "\"trap_masked\":%s,\"errno_seed\":%d,\"ambient_restored\":%s,\"setup_error\":%s,"
                "\"first_mismatches\":[", passed ? "true" : "false", input_count,
                static_cast<unsigned long long>(comparisons), static_cast<unsigned long long>(comparisons * 2),
                static_cast<unsigned long long>(mismatches), mxcsr_mask, trap_masked ? "true" : "false",
                EDOM, restored ? "true" : "false",
                setup_error ? "\"setup_or_probe_exception\"" : "null");
    for (std::size_t i = 0; i < failures.size(); ++i) {
        const auto& failure = failures[i];
        std::printf("%s{\"input\":\"%08x\",\"round\":%d,\"FTZ_DAZ\":%u,\"preexisting_inexact\":%d,"
                    "\"reference\":", i ? "," : "", failure.input, failure.rounding, failure.controls, failure.seed);
        print_observed(failure.reference); std::printf(",\"candidate\":");
        print_observed(failure.candidate); std::printf("}");
    }
    std::printf("]}\n");
    const bool output_ok = !std::ferror(stdout);
    const bool final_restore = restore();
    return passed && output_ok && final_restore ? 0 : 1;
}
'''

if __name__ == '__main__':
    proof = json.loads((W.parent / 'proof.json').read_text())
    assert pin(SOURCE) == {k: proof['source'][k] for k in ('bytes', 'sha256')}
    assert pin(PATCH) == {k: proof['host_patch'][k] for k in ('bytes', 'sha256')}
    old = SOURCE.read_text(); new = apply(old, PATCH.read_text())
    reference, candidate = spacing(old), spacing(new)
    assert old.replace(reference, candidate, 1) == new
    assert reference != candidate and CPP.count('@REFERENCE@') == CPP.count('@CANDIDATE@') == 1
    generated = CPP.replace('@REFERENCE@', reference).replace('@CANDIDATE@', candidate)
    with (W / 'probe.cpp').open('x') as stream:
        stream.write(generated)
    receipt = {'source': str(SOURCE.relative_to(R)), 'source_identity': pin(SOURCE),
               'patch': str(PATCH.relative_to(R)), 'patch_identity': pin(PATCH),
               'reference_lambda_sha256': hashlib.sha256(reference.encode()).hexdigest(),
               'candidate_lambda_sha256': hashlib.sha256(candidate.encode()).hexdigest(),
               'generated_source': pin(W / 'probe.cpp'), 'generator': pin(Path(__file__)),
               'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
               'scope': 'Actual source/patched lambdas; no formula substitute, no compilation or execution.'}
    with (W / 'generation.json').open('x') as stream:
        json.dump(receipt, stream, indent=2); stream.write('\n')
