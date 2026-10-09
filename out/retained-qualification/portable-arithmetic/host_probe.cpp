#define PORTABLE_ARITHMETIC_HOST 1
#include "portable_f32.h"
#include <iostream>
#include <iomanip>
static_assert(sizeof(uint) == 4 && sizeof(int) == 4);
int main() {
    uint operation, a, b;
    while (std::cin >> std::hex >> operation >> a >> b)
        std::cout << std::hex << std::setfill('0') << std::setw(8)
                  << PaEvaluate(operation, a, b) << '\n';
    return std::cin.eof() ? 0 : 1;
}
