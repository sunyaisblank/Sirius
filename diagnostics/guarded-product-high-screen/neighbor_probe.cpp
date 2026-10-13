#include "portable_binary32.h"
#include <fstream>
#include <iostream>
#include <array>

namespace sirius::portable_binary32 {
uint PB32HighModel = 0;
uint PB32GuardedCalls = 0;
uint PB32Corrections = 0;
uint PB32DownwardBinades = 0;
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    std::ifstream input(argv[1], std::ios::binary);
    if (!input) return 3;
    using namespace sirius::portable_binary32;
    std::array<uint, 7> row{};
    std::size_t count = 0;
    while (input.read(reinterpret_cast<char*>(row.data()), sizeof(row))) {
        PB32HighModel = row[2];
        uint before = PB32GuardedCalls;
        auto result = PB32ProductResidual(row[0], row[1]);
        if (result.high != row[3] || result.low != row[4] || result.valid != row[5] ||
            PB32GuardedCalls - before != row[6]) {
            std::cerr << "mismatch record " << count << ' ' << std::hex
                      << row[0] << ' ' << row[1] << ' ' << row[2] << ' '
                      << result.high << ' ' << result.low << ' ' << result.valid << '\n';
            return 4;
        }
        ++count;
    }
    if (!input.eof() || input.gcount() != 0) return 5;
    std::cout << "{\"records\":" << count << ",\"guarded_calls\":" << PB32GuardedCalls
              << ",\"corrections\":" << PB32Corrections
              << ",\"downward_binades\":" << PB32DownwardBinades << "}\n";
}
