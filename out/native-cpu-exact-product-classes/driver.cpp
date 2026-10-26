// Fixed mathematical centers plus exact complete-bank comparison with sealed v3.
// No candidate predicate, random generator, scalar FP oracle or benchmark.
#include <array>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <immintrin.h>

extern "C" void sirius_native_dense_multiply8(std::uint32_t*, std::uint32_t,
    std::uint32_t, std::uint32_t, std::uint32_t) noexcept;
extern "C" void sirius_native_v3_dense_multiply8(std::uint32_t*, std::uint32_t,
    std::uint32_t, std::uint32_t, std::uint32_t) noexcept;

namespace {
struct Lane { std::uint32_t a, b, high, low; };
struct Case {
    const char* name;
    bool independent_center;
    std::array<Lane, 8> lanes;
};

constexpr Case Uniform(const char* name, std::uint32_t a, std::uint32_t b,
                       bool independent = false, std::uint32_t high = 0,
                       std::uint32_t low = 0) {
    Case result{name, independent, {}};
    result.lanes.fill({a, b, high, low});
    return result;
}

constexpr auto cases = [] {
    std::array<Case, 17> result{{
        Uniform("A_coarse_positive", 0x3fc00000U, 0x3fe00000U, true, 0x40280000U),
        Uniform("A_coarse_all_signs", 0x3fc00000U, 0x3fe00000U, true, 0x40280000U),
        Uniform("A_joint_floor", 0x20000000U, 0x20000000U, true, 0x00800000U),
        Uniform("A_below_joint_floor_reference_only", 0x1f800000U, 0x20000000U),
        Uniform("B_mixed_power_side_and_all_signs", 0x3f800fffU, 0x40000000U, true, 0x40000fffU),
        Uniform("B_joint_floor", 0x4b000fffU, 0x00800000U, true, 0x0c000fffU),
        Uniform("B_below_joint_floor_reference_only", 0x4a800fffU, 0x00800000U),
        Uniform("B_other_exponent_floor", 0x0c000fffU, 0x5d800000U, true, 0x2a000fffU),
        Uniform("B_below_other_floor_reference_only", 0x0b800fffU, 0x5d800000U),
        Uniform("B_subnormal_power_reference_only", 0x00000001U, 0x56800fffU),
        Uniform("Double_one_ulp", 0x3f800001U, 0x3f800001U, true, 0x3f800002U, 0x28800000U),
        Uniform("Double_joint_floor", 0x2b800001U, 0x2b800001U, true, 0x17800002U, 0x00800000U),
        Uniform("Double_below_joint_floor_reference_only", 0x2b000001U, 0x2b800001U),
        Uniform("exact_nonpower_zero_residual", 0x3fc00000U, 0x3f800ffeU, true, 0x3fc017fdU),
        Uniform("maximum_operands_reference_only", 0x5d800000U, 0x5d800000U),
        Uniform("above_operand_bound_reference_only", 0x5d800001U, 0x5d800000U),
        Uniform("zero_operands_reference_only", 0, 0x3f800fffU),
    }};
    for (unsigned lane = 0; lane < 8; ++lane) {
        const auto sign_a = (lane & 1U) << 31;
        const auto sign_b = ((lane >> 1) & 1U) << 31;
        result[1].lanes[lane] = {0x3fc00000U ^ sign_a, 0x3fe00000U ^ sign_b,
                                 0x40280000U ^ sign_a ^ sign_b, 0};
        const auto power_sign = ((lane >> 1) & 1U) << 31;
        const auto other_sign = ((lane >> 2) & 1U) << 31;
        const auto power = 0x40000000U ^ power_sign;
        const auto other = 0x3f800fffU ^ other_sign;
        result[4].lanes[lane] = {(lane & 1U) ? power : other,
                                 (lane & 1U) ? other : power,
                                 0x40000fffU ^ power_sign ^ other_sign, 0};
        result[16].lanes[lane] = {(lane & 1U) ? 0x80000000U : 0,
                                  (lane & 2U) ? 0xbf800fffU : 0x3f800fffU, 0, 0};
    }
    return result;
}();

// Exact integer dyadics verify the fixed center constants. This performs no
// FP arithmetic and does not classify operands using any candidate guard.
struct Dyadic { std::int64_t integer; int exponent; };
constexpr Dyadic Decode(std::uint32_t bits) {
    const auto biased = (bits >> 23) & 255U;
    std::uint32_t magnitude = (bits & 0x007fffffU) | (biased ? 0x00800000U : 0);
    int exponent = biased ? static_cast<int>(biased) - 150 : -149;
    while (magnitude && !(magnitude & 1U)) { magnitude >>= 1; ++exponent; }
    const auto integer = static_cast<std::int64_t>(magnitude);
    return {(bits >> 31) ? -integer : integer, exponent};
}
constexpr bool Scale(Dyadic value, int exponent, std::int64_t& result) {
    if (!value.integer) { result = 0; return true; }
    const int shift = value.exponent - exponent;
    if (shift < 0 || shift > 62) return false;
    const auto magnitude = static_cast<std::uint64_t>(value.integer < 0 ? -value.integer : value.integer);
    if (magnitude > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()) >> shift)
        return false;
    const auto scaled = static_cast<std::int64_t>(magnitude << shift);
    result = value.integer < 0 ? -scaled : scaled;
    return true;
}
constexpr bool ExactCenter(const Lane& lane) {
    for (auto bits : {lane.a, lane.b, lane.high, lane.low})
        if (((bits >> 23) & 255U) == 255U) return false;
    const auto a = Decode(lane.a), b = Decode(lane.b);
    const Dyadic product{a.integer * b.integer, a.exponent + b.exponent};
    const auto high = Decode(lane.high), low = Decode(lane.low);
    int exponent = product.exponent;
    if (high.integer && high.exponent < exponent) exponent = high.exponent;
    if (low.integer && low.exponent < exponent) exponent = low.exponent;
    std::int64_t p = 0, h = 0, l = 0;
    return Scale(product, exponent, p) && Scale(high, exponent, h) && Scale(low, exponent, l)
        && p == h + l;
}
constexpr bool AllCentersExact() {
    unsigned independent = 0;
    for (const auto& test : cases) if (test.independent_center) {
        ++independent;
        for (const auto& lane : test.lanes) if (!ExactCenter(lane)) return false;
    }
    return independent == 9;
}
static_assert(AllCentersExact());

bool CenterWord(std::uint32_t actual, std::uint32_t expected) {
    // Independent real-number equality ignores zero sign; complete-bank
    // comparison below remains bit-exact for every sign and every field.
    return expected == 0 ? (actual & 0x7fffffffU) == 0 : actual == expected;
}
} // namespace

int main() {
    const bool v3 = __builtin_cpu_supports("x86-64-v3");
    const bool v4 = __builtin_cpu_supports("x86-64-v4");
    if (!v3 || !v4) {
        std::cout << "{\"passed\":false,\"unsupported_actual_CPU\":true}\n";
        return 2;
    }
    constexpr unsigned prefix = 8, bank_words = 120, suffix = 8;
    using Storage = std::array<std::uint32_t, prefix + bank_words + suffix>;
    const unsigned saved = _mm_getcsr();
    const unsigned scope = (saved | 0x1f80U) & ~(0xe040U | 0x3fU);
    bool passed = true, all_equal = true, all_immutable = true, all_canaries = true;
    bool all_centers = true, all_controls = true, all_sticky = true, all_restored = true;
    bool subnormal_reference_denormal_operand = false;
    std::cout << std::boolalpha << "{\"cases\":[";
    for (unsigned index = 0; index < cases.size(); ++index) {
        const auto& test = cases[index];
        Storage actual;
        actual.fill(0xa5a5a5a5U);
        auto* bank = actual.data() + prefix;
        for (unsigned lane = 0; lane < 8; ++lane) {
            bank[lane] = test.lanes[lane].a;
            bank[40 + lane] = test.lanes[lane].b;
            for (unsigned reg = 0; reg < 2; ++reg) {
                for (unsigned field = 1; field < 4; ++field) bank[reg * 40 + field * 8 + lane] = 0;
                bank[reg * 40 + 32 + lane] = 1;
            }
        }
        const Storage input = actual;
        Storage reference = actual;
        _mm_setcsr(scope);
        sirius_native_dense_multiply8(bank, 2, 0, 1, 255);
        const unsigned candidate_mxcsr = _mm_getcsr();
        _mm_setcsr(scope);
        sirius_native_v3_dense_multiply8(reference.data() + prefix, 2, 0, 1, 255);
        const unsigned reference_mxcsr = _mm_getcsr();
        const bool controls = (candidate_mxcsr & ~0x3fU) == scope && (reference_mxcsr & ~0x3fU) == scope;
        const bool sticky = candidate_mxcsr == reference_mxcsr;
        _mm_setcsr(saved);
        const bool restored = _mm_getcsr() == saved;
        bool equal = true, immutable = true, canaries = true, center = true;
        for (unsigned i = 0; i < bank_words; ++i)
            equal &= actual[prefix + i] == reference[prefix + i];
        for (unsigned i = 0; i < 80; ++i)
            immutable = immutable && bank[i] == input[prefix + i] && reference[prefix + i] == input[prefix + i];
        for (unsigned i = 0; i < prefix; ++i)
            canaries = canaries && actual[i] == input[i] && reference[i] == input[i];
        for (unsigned i = prefix + bank_words; i < actual.size(); ++i)
            canaries = canaries && actual[i] == input[i] && reference[i] == input[i];
        if (test.independent_center) for (unsigned lane = 0; lane < 8; ++lane)
            center = center && CenterWord(bank[80 + lane], test.lanes[lane].high)
                && CenterWord(bank[88 + lane], test.lanes[lane].low)
                && CenterWord(bank[96 + lane], 0) && bank[112 + lane] == 1;
        if (index == 9) {
            subnormal_reference_denormal_operand = (reference_mxcsr & 2U) != 0;
            passed = passed && subnormal_reference_denormal_operand;
        }
        const bool case_passed = equal && immutable && canaries && center && controls && sticky && restored
            && (index != 9 || subnormal_reference_denormal_operand);
        passed = passed && case_passed;
        all_equal = all_equal && equal; all_immutable = all_immutable && immutable;
        all_canaries = all_canaries && canaries; all_centers = all_centers && center;
        all_controls = all_controls && controls; all_sticky = all_sticky && sticky;
        all_restored = all_restored && restored;
        if (index) std::cout << ',';
        std::cout << "{\"name\":\"" << test.name << "\",\"passed\":" << case_passed
            << ",\"independent_center_checked\":" << test.independent_center
            << ",\"whole_bank_equal\":" << equal << ",\"input_immutable\":" << immutable
            << ",\"canaries_preserved\":" << canaries << ",\"independent_center_equal\":" << center
            << ",\"MXCSR_controls_unchanged\":" << controls << ",\"MXCSR_sticky_equal\":" << sticky
            << ",\"candidate_MXCSR\":" << candidate_mxcsr << ",\"reference_MXCSR\":" << reference_mxcsr
            << ",\"caller_MXCSR_state_restored\":" << restored << '}';
    }
    std::cout << "],\"passed\":" << passed
        << ",\"case_count\":17,\"candidate_calls\":17,\"original_v3_reference_calls\":17"
        << ",\"active_lane_results\":136,\"compared_bank_words\":2040,\"compared_canary_words\":272"
        << ",\"immutable_input_words_per_case\":80,\"independent_center_cases\":9"
        << ",\"independent_center_lane_results\":72,\"reference_only_cases\":8"
        << ",\"whole_bank_equal\":" << all_equal << ",\"input_immutable\":" << all_immutable
        << ",\"canaries_preserved\":" << all_canaries << ",\"independent_center_equal\":" << all_centers
        << ",\"MXCSR_controls_unchanged\":" << all_controls << ",\"MXCSR_sticky_equal\":" << all_sticky
        << ",\"caller_MXCSR_state_restored\":" << all_restored
        << ",\"subnormal_power_reference_denormal_operand_observed\":" << subnormal_reference_denormal_operand
        << ",\"actual_CPU_v3\":" << v3 << ",\"actual_CPU_v4\":" << v4
        << ",\"scope\":\"Fixed 17-packet actual-dispatcher/sealed-v3 comparison; nine integer-dyadic center identities only. Other cases assert reference equivalence only. No random cases, benchmark, branch-PC coverage, hit-rate, provider/frame or scientific/release acceptance.\"}\n";
    return passed ? 0 : 3;
}
