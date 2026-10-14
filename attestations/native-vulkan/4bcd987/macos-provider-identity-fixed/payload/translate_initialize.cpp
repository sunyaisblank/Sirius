// Finite issue-48 CPU translation only. No Vulkan, Metal, compiler subprocess or dispatch.
#include "spirv_msl.hpp"
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::runtime_error("usage: translate_initialize input.spv fresh-output.metal");
        const std::filesystem::path input(argv[1]), output(argv[2]);
        if (std::filesystem::exists(output) || std::filesystem::is_symlink(output))
            throw std::runtime_error("output must be fresh");
        std::ifstream stream(input, std::ios::binary);
        if (!stream) throw std::runtime_error("input could not be opened");
        const std::vector<char> bytes((std::istreambuf_iterator<char>(stream)), {});
        if (bytes.size() != 77444) throw std::runtime_error("different captured SPIR-V size");
        std::vector<std::uint32_t> words(bytes.size() / 4);
        std::memcpy(words.data(), bytes.data(), bytes.size());
        if (words[0] != 0x07230203u) throw std::runtime_error("different SPIR-V byte order or magic");
        spirv_cross::CompilerMSL compiler(std::move(words));
        compiler.set_entry_point("main", spv::ExecutionModelGLCompute);
        auto options = compiler.get_msl_options();
        options.platform = spirv_cross::CompilerMSL::Options::macOS;
        options.msl_version = 30200;
        options.argument_buffers = true;
        options.argument_buffers_tier = spirv_cross::CompilerMSL::Options::ArgumentBuffersTier::Tier2;
        options.pad_argument_buffer_resources = true;
        options.force_active_argument_buffer_resources = false;
        options.buffer_size_buffer_index = 29;
        options.dynamic_offsets_buffer_index = 30;
        options.indirect_params_buffer_index = 28;
        options.texture_buffer_native = true;
        options.dispatch_base = false;
        options.replace_recursive_inputs = true;
        options.fixed_subgroup_size = 32;
        options.emulate_subgroups = false;
        compiler.set_msl_options(options);
        auto common = compiler.get_common_options();
        common.vertex.flip_vert_y = true;
        common.vertex.fixup_clipspace = true;
        compiler.set_common_options(common);
        for (std::uint32_t index = 0; index < 3; ++index) {
            spirv_cross::MSLResourceBinding binding;
            binding.stage = spv::ExecutionModelGLCompute;
            binding.basetype = spirv_cross::SPIRType::Void;
            binding.desc_set = 0;
            binding.binding = index == 0 ? spirv_cross::kBufferSizeBufferBinding : index - 1;
            binding.count = 1;
            binding.msl_buffer = index;
            compiler.add_msl_resource_binding(binding);
        }
        const auto msl = compiler.compile();
        std::ofstream result(output, std::ios::binary | std::ios::out);
        if (!result) throw std::runtime_error("output could not be opened");
        result.write(msl.data(), static_cast<std::streamsize>(msl.size()));
        result.close();
        if (!result) throw std::runtime_error("short output write");
        std::cerr << "CPU translation completed using a source-backed configuration reconstruction; exact byte acceptance remains separate.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
