// Contingent packing adapter only. Dense arithmetic is reused byte-for-byte.
// Inputs and output are initialized five-column, eight-lane word buffers.
// Snapshot all three before writes, so complete or partial aliases are safe.
#include <cstdint>
#include <cstring>

extern "C" void sirius_native_dense_add8(std::uint32_t*,std::uint32_t,std::uint32_t,std::uint32_t,std::uint32_t) noexcept;
extern "C" void sirius_native_dense_multiply8(std::uint32_t*,std::uint32_t,std::uint32_t,std::uint32_t,std::uint32_t) noexcept;

extern "C" __attribute__((visibility("hidden")))
void sirius_native_dense_values_add8(const std::uint32_t* a,const std::uint32_t* b,
                                    std::uint32_t* output,std::uint32_t mask) noexcept {
    alignas(32) std::uint32_t bank[120];
    std::memcpy(bank,a,160);
    std::memcpy(bank+40,b,160);
    std::memcpy(bank+80,output,160);
    sirius_native_dense_add8(bank,2,0,1,mask);
    std::memcpy(output,bank+80,160);
}

extern "C" __attribute__((visibility("hidden")))
void sirius_native_dense_values_multiply8(const std::uint32_t* a,const std::uint32_t* b,
                                         std::uint32_t* output,std::uint32_t mask) noexcept {
    alignas(32) std::uint32_t bank[120];
    std::memcpy(bank,a,160);
    std::memcpy(bank+40,b,160);
    std::memcpy(bank+80,output,160);
    sirius_native_dense_multiply8(bank,2,0,1,mask);
    std::memcpy(output,bank+80,160);
}
