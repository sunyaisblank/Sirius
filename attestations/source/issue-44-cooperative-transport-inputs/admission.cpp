#include "sirius/backend/device.h"
#include "retained_kernels.h"
#include "original_transport.h"
#include "support/retained_transport/reference_cases.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

using std::uint32_t;
using namespace sirius::backend;
using namespace sirius::backend::retained_program;
constexpr std::size_t rows = 8;
using Packet = std::array<uint32_t, 230>;

void require(bool value, const std::string& message) {
    if (!value) throw std::runtime_error(message);
}

void save(const std::filesystem::path& path, std::span<const uint32_t> words) {
    std::ofstream stream(path, std::ios::binary);
    stream.write(reinterpret_cast<const char*>(words.data()),
                 static_cast<std::streamsize>(words.size_bytes()));
    require(bool(stream), "write capture " + path.string());
}

int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: admission owned-output-folder");
        const auto folder = std::filesystem::path(argv[1]) / "captures";
        std::filesystem::create_directory(folder);
        const auto inventory = EnumerateVulkanDevices();
        require(bool(inventory), "enumerate software provider");
        const auto selected = ResolveVulkanDeviceIndex(*inventory);
        require(bool(selected), "select software provider");
        namespace fixture = sirius::test::retained_transport;
        require(std::string(fixture::cases.front().name) == "critical-0" &&
                    std::string(fixture::cases.back().name) == "analytic-flat",
                "unchanged curved/flat fixtures");
        const std::array good{fixture::cases.front().input, fixture::cases.back().input};
        std::vector<Packet> malformed;
        for (std::size_t field = 0; field < 46; ++field) {
            auto invalid = good[0];
            invalid[5 * field + 4] = 0;
            malformed.push_back(invalid);
            invalid = good[0];
            invalid[5 * field + 3] = 0xbf800000U;
            malformed.push_back(invalid);
        }
        for (const std::size_t field : {0U, 15U, 16U, 31U, 32U, 44U, 45U}) {
            for (unsigned kind = 0; kind < 5; ++kind) {
                auto invalid = good[0];
                const std::array<unsigned, 5> component{4, 0, 1, 2, 0};
                const std::array<uint32_t, 5> word{2, 0x7fc00000U, 0x7f800000U,
                                                 0x7fc00000U, 0x7c000000U};
                invalid[5 * field + component[kind]] = word[kind];
                malformed.push_back(invalid);
            }
        }
        require(malformed.size() == 127, "complete malformed partitions");
        const std::array<std::span<const uint32_t>, 3> original{
            diagnostic_original::kTransportPortableShader,
            diagnostic_original::kTransportPortableFp64Shader,
            diagnostic_original::kTransportPortableNormalSumShader};
        const std::array<std::span<const uint32_t>, 3> candidate{
            kTransportPortableShader, kTransportPortableFp64Shader,
            kTransportPortableNormalSumShader};
        std::size_t dispatches = 0, invalid_observations = 0, compared_words = 0;
        for (std::size_t mode = 0; mode < candidate.size(); ++mode) {
            auto opened = CreateVulkanDevice(*selected);
            require(bool(opened), "open actual software device");
            auto& device = **opened;
            const auto& info = device.Info();
            require(info.driver_id == 13 && info.driver_name == "llvmpipe" &&
                        info.kind == DeviceKind::kSoftware && info.rounds_fp32_to_nearest &&
                        info.supports_fp64 && info.rounds_fp64_to_nearest,
                    "actual pinned software/provider/product modes");
            require(bool(device.SetBufferAllocationLimit(8 * 1024 * 1024)), "8 MiB limit");
            const auto input = device.CreateBuffer((1 + 230 * rows + kTransportProgram.size()) * 4,
                                                    BufferUsage::kStorage);
            const auto output = device.CreateBuffer(rows * kTransportRowWords * 4,
                                                     BufferUsage::kStorage);
            require(bool(input) && bool(output), "original transport buffers");
            const std::array buffers{*input, *output};
            std::array<KernelHandle, 2> kernels{};
            for (std::size_t producer = 0; producer < 2; ++producer) {
                const auto code = producer == 0 ? original[mode] : candidate[mode];
                const auto loaded = device.LoadKernel(code);
                require(bool(loaded), "load complete original/candidate module");
                kernels[producer] = *loaded;
                save(folder / (std::to_string(mode) + "-module-" +
                               std::to_string(producer) + ".spv"), code);
            }
            const auto run = [&](const std::array<Packet, rows>& packets, std::size_t number) {
                std::vector<uint32_t> words(1 + 230 * rows + kTransportProgram.size());
                words[0] = rows;
                for (std::size_t row = 0; row < rows; ++row)
                    std::copy(packets[row].begin(), packets[row].end(), words.begin() + 1 + 230 * row);
                std::copy(kTransportProgram.begin(), kTransportProgram.end(),
                          words.begin() + 1 + 230 * rows);
                const auto prefix = std::to_string(mode) + "-" + std::to_string(number);
                save(folder / (prefix + "-input.bin"), words);
                std::array<std::vector<uint32_t>, 2> result;
                for (std::size_t producer = 0; producer < 2; ++producer) {
                    require(bool(device.WriteBuffer(*input, std::as_bytes(std::span(words)))),
                            "write original complete input");
                    result[producer].assign(rows * kTransportRowWords, 0xdeadbeefU);
                    require(bool(device.WriteBuffer(*output, std::as_bytes(std::span(result[producer])))),
                            "poison complete output");
                    require(bool(device.Dispatch(kernels[producer], buffers, rows, 1, 1)),
                            "complete transport dispatch");
                    ++dispatches;
                    require(bool(device.ReadBuffer(*output,
                                                   std::as_writable_bytes(std::span(result[producer])))),
                            "complete transport readback");
                    save(folder / (prefix + "-output-" + std::to_string(producer) + ".bin"),
                         result[producer]);
                }
                require(result[0] == result[1], "complete same-provider original/candidate words");
                compared_words += result[0].size();
                return result[1];
            };
            std::array<Packet, rows> packets{};
            for (std::size_t row = 0; row < rows; ++row) packets[row] = good[row % 2];
            const auto reference = run(packets, 0);
            for (std::size_t row = 0; row < rows; ++row)
                require(reference[row * kTransportRowWords] == 1 &&
                            reference[row * kTransportRowWords + 1] == 7,
                        "original valid curved/flat seven-stage reference");
            std::size_t number = 1;
            for (std::size_t start = 0; start < malformed.size(); start += rows - 2) {
                packets.fill(good[0]);
                packets.back() = good[1];
                const auto count = std::min(rows - 2, malformed.size() - start);
                for (std::size_t i = 0; i < count; ++i) packets[1 + i] = malformed[start + i];
                const auto actual = run(packets, number++);
                for (std::size_t row = 0; row < rows; ++row) {
                    const auto first = actual.begin() + row * kTransportRowWords;
                    const auto last = first + kTransportRowWords;
                    if (row >= 1 && row <= count) {
                        require(std::all_of(first, last, [](uint32_t word) { return word == 0; }),
                                "every malformed field preserves the complete zero row");
                        invalid_observations += 2;
                    } else {
                        const auto expected = reference.begin() + (row == rows - 1 ? 1 : 0) * kTransportRowWords;
                        require(std::equal(first, last, expected), "valid neighboring row isolation");
                    }
                }
            }
            for (std::size_t row = 0; row < rows; ++row) packets[row] = good[row % 2];
            require(run(packets, number) == reference, "complete valid recovery after all refusals");
            require(device.BufferAllocationBytes() <= 8 * 1024 * 1024,
                    "preserved original device allocation limit");
        }
        std::cout << "{\"passed\":true,\"scope\":\"field admission, complete original words and recovery; no independent full science or benefit claim\","
                  << "\"product_scalar_variants\":3,\"malformed_partitions_per_variant\":127,"
                  << "\"invalid_row_observations\":" << invalid_observations
                  << ",\"dispatches\":" << dispatches << ",\"complete_compared_words\":"
                  << compared_words << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
