#include "portable_binary32.h"
#include "word_io.h"
#include <array>
#include <iostream>
#include <limits>

using namespace sirius::portable_binary32;

std::uint64_t IntegerRoot(std::uint64_t n) {
    std::uint64_t low = 0, high = std::uint64_t(1) << 28;
    while (low + 1 < high) {
        const auto mid = (low + high) / 2;
        if (mid * mid <= n) low = mid;
        else high = mid;
    }
    return low;
}

std::uint32_t Evaluate(std::uint32_t op, std::uint32_t a, std::uint32_t b) {
    switch (op) {
        case 0: return PB32Add(a,b);
        case 1: return PB32Subtract(a,b);
        case 2: return PB32Multiply(a,b);
        case 3: return PB32Divide(a,b);
        case 4: return PB32Sqrt(a);
        case 5: return PB32Equal(a,b) ? 1u : 0u;
        case 6: return PB32Less(a,b) ? 1u : 0u;
        case 7: return PB32LessEqual(a,b) ? 1u : 0u;
        case 8: return PB32FromUnsigned(a);
        case 9: return PB32FromSignedBits(a);
        case 10: return PB32Abs(a);
        case 11: return PB32Negate(a);
        case 12: return PB32Classify(a);
        default: throw std::runtime_error("invalid operation");
    }
}

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::runtime_error("input expected required");
        const auto input = ReadWords(argv[1]), expected = ReadWords(argv[2]);
        CheckShape(input, expected);
        std::uint64_t accepted = 0, refused = 0, checked = 0;
        std::uint32_t random = 0x7249a1e3u;
        auto next = [&] { random ^= random << 13; random ^= random >> 17; random ^= random << 5; return random; };
        for (unsigned i = 0; i < 20008; ++i) {
            const std::array boundary{0x800000u,0x800001u,0xffffffu,0xfffffeu};
            const auto d = i < 8 ? boundary[i % 4] : 0x800000u | (next() & 0x7fffffu);
            const auto r = i < 4 ? d : i < 8 ? 2*d-1 : d + next() % d;
            const auto n = std::uint64_t(r) << 26;
            const auto shift = 29u + (i & 1u);
            const auto sn = std::uint64_t(d) << shift;
            PB32Wide wide{std::uint32_t(sn),std::uint32_t(sn >> 32)};
            const auto dq = n/d, sq = IntegerRoot(sn);
            for (unsigned kind = 0; kind < 2; ++kind) {
                const auto q = kind == 0 ? dq : sq;
                const std::array seeds{double(q)-3.75,double(q)-2.0,double(q)-1.0,double(q)-0.75,
                    double(q),double(q)+0.75,double(q)+1.0,double(q)+1.75,double(q)+2.0,
                    double(q)+100.0,0.0,-0.0,-1.0,268435456.0,
                    std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::denorm_min()};
                for (double seed : seeds) {
                    const auto result = kind == 0 ? PB32CheckedDivide(r,d,seed) : PB32CheckedSqrt(wide,seed);
                    ++checked;
                    if (result.valid == 0) { ++refused; continue; }
                    ++accepted;
                    const auto residue = kind == 0 ? n - q*d : sn-q*q;
                    if (result.quotient != q || result.sticky != (residue != 0 ? 1u : 0u))
                        throw std::runtime_error("exact postcondition admitted wrong quotient/sticky");
                }
            }
        }
        if (accepted == 0 || refused == 0) throw std::runtime_error("vacuous seed partition");
        for (std::size_t i = 0; i < input[0]; ++i) {
            const auto op=input[1+3*i], a=input[2+3*i], b=input[3+3*i];
            PB32Product result = op == 13 ? PB32ProductResidual(a,b) : PB32Product{Evaluate(op,a,b),0,1};
            if (result.high != expected[3*i] || result.low != expected[3*i+1] || result.valid != expected[3*i+2])
                throw std::runtime_error("independent original corpus mismatch at " + std::to_string(i));
        }
        std::cout << "{\"passed\":true,\"injected_seeds\":" << checked
            << ",\"accepted_exact\":" << accepted << ",\"refused\":" << refused
            << ",\"original_records\":" << input[0] << ",\"corpus_mismatches\":0}\n";
        return 0;
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
