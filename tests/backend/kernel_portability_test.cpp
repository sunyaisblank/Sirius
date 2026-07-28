// Portability gate (docs/ARCHITECTURE.md section 4): the single Slang trace
// source emits vendor-native kernel source alongside its SPIR-V artefacts —
// CUDA C++ for the nvcc toolchain and Metal Shading Language for Apple
// silicon. Neither toolchain exists on the measured machines, so what this
// suite pins is the emission contract: the artefact exists, is substantial,
// and carries the compute entry point in the vendor's native spelling. A
// kernel edit that breaks a vendor target fails this gate rather than
// surfacing on first contact with that vendor's hardware.

#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <string>

namespace {

[[nodiscard]] std::string ReadArtefact(const std::string& name) {
#ifndef SIRIUS_KERNEL_DIR
    return {};
#else
    std::ifstream file(std::string(SIRIUS_KERNEL_DIR) + "/portability/" + name);
    if (!file) {
        return {};
    }
    std::ostringstream contents;
    contents << file.rdbuf();
    return contents.str();
#endif
}

TEST(KernelPortability, CudaEmissionCarriesTheNativeComputeEntryPoint) {
#ifndef SIRIUS_KERNEL_DIR
    GTEST_SKIP() << "kernels not compiled (slangc absent at configure time)";
#endif
    const std::string cuda = ReadArtefact("trace.cu");
    ASSERT_FALSE(cuda.empty()) << "trace.cu missing or unreadable";
    EXPECT_NE(cuda.find("__global__ void ComputeMain"), std::string::npos)
        << "CUDA emission lost its __global__ compute entry point";
}

TEST(KernelPortability, MetalEmissionCarriesTheNativeComputeEntryPoint) {
#ifndef SIRIUS_KERNEL_DIR
    GTEST_SKIP() << "kernels not compiled (slangc absent at configure time)";
#endif
    const std::string metal = ReadArtefact("trace.metal");
    ASSERT_FALSE(metal.empty()) << "trace.metal missing or unreadable";
    EXPECT_NE(metal.find("[[kernel]] void ComputeMain"), std::string::npos)
        << "Metal emission lost its [[kernel]] compute entry point";
    EXPECT_NE(metal.find("thread_position_in_grid"), std::string::npos)
        << "Metal emission lost its grid-position binding";
}

}  // namespace
