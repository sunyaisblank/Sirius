#include "fp64_scalar.h"
#include "word_io.h"
#include <array>
#include <cfenv>
#include <iostream>
int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::runtime_error("usage: host_probe input expected actual");
        if (std::fegetround() != FE_TONEAREST) throw std::runtime_error("host RTE required");
        const auto input = ReadWords(argv[1]), expected = ReadWords(argv[2]);
        CheckShape(input, expected);
        std::vector<std::uint32_t> actual(expected.size());
        std::size_t mismatches = 0;
        for (std::size_t row = 0; row < input[0]; ++row) {
            const auto result = sirius::portable_fp64_scalar::PB64Evaluate(
                input[1 + 3 * row], input[2 + 3 * row], input[3 + 3 * row]);
            const std::array words{result.high, result.low, result.valid};
            for (std::size_t column = 0; column < 3; ++column) {
                actual[3 * row + column] = words[column];
                if (words[column] != expected[3 * row + column]) {
                    if (mismatches < 8) std::cerr << "mismatch row=" << row << " col=" << column << '\n';
                    ++mismatches;
                }
            }
        }
        WriteWords(argv[3], actual);
        std::cout << "{\"scope\":\"isolated host raw-word prototype\",\"cases\":" << input[0]
                  << ",\"observed_words\":" << actual.size() << ",\"mismatches\":" << mismatches << "}\n";
        return mismatches == 0 ? 0 : 1;
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
