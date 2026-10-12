#pragma once
#include <cstdint>
#include <fstream>
#include <span>
#include <stdexcept>
#include <vector>
inline std::vector<std::uint32_t> ReadWords(const char* path) {
    std::ifstream stream(path, std::ios::binary | std::ios::ate);
    const auto size = stream.tellg();
    if (!stream || size <= 0 || size > 8 * 1024 * 1024 || size % 4 != 0)
        throw std::runtime_error("invalid/bounded word input");
    std::vector<unsigned char> bytes(static_cast<std::size_t>(size));
    stream.seekg(0);
    if (!stream.read(reinterpret_cast<char*>(bytes.data()), size))
        throw std::runtime_error("input read failed");
    std::vector<std::uint32_t> words(bytes.size() / 4);
    for (std::size_t i = 0; i < words.size(); ++i)
        for (unsigned j = 0; j < 4; ++j)
            words[i] |= std::uint32_t(bytes[4 * i + j]) << (8 * j);
    return words;
}
inline void WriteWords(const char* path, std::span<const std::uint32_t> words) {
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("output open failed");
    for (auto word : words)
        for (unsigned j = 0; j < 4; ++j) stream.put(char((word >> (8 * j)) & 255u));
    stream.close();
    if (!stream) throw std::runtime_error("output close failed");
}
inline void CheckShape(const std::vector<std::uint32_t>& input,
                       const std::vector<std::uint32_t>& expected) {
    if (input.empty() || input[0] == 0u || input[0] > 262143u ||
        input.size() != 1 + 4 * std::size_t(input[0]) ||
        expected.size() != 4 * std::size_t(input[0]))
        throw std::runtime_error("bounded four-word scalar corpus shape mismatch");
}
