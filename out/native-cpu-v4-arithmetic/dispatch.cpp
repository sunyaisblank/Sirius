// Compile without optional ISA flags. The provider retains its original v3 guard.
// No vector or Triple crosses this original pointer/u32 ABI.
#include <cstdint>
extern "C" __attribute__((visibility("hidden")))
bool sirius_native_v4_supported() noexcept {
    return __builtin_cpu_supports("x86-64-v4");
}
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v3_dense_add8(std::uint32_t* scratch,std::uint32_t destination,
                             std::uint32_t a_index,std::uint32_t b_index,
                             std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v4_dense_add8(std::uint32_t* scratch,std::uint32_t destination,
                             std::uint32_t a_index,std::uint32_t b_index,
                             std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_dense_add8(std::uint32_t* scratch,std::uint32_t destination,
                             std::uint32_t a_index,std::uint32_t b_index,
                             std::uint32_t mask) noexcept {
    if(sirius_native_v4_supported()) sirius_native_v4_dense_add8(scratch,destination,a_index,b_index,mask);
    else sirius_native_v3_dense_add8(scratch,destination,a_index,b_index,mask);
}
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v3_dense_multiply8(std::uint32_t* scratch,std::uint32_t destination,
                                  std::uint32_t a_index,std::uint32_t b_index,
                                  std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v4_dense_multiply8(std::uint32_t* scratch,std::uint32_t destination,
                                  std::uint32_t a_index,std::uint32_t b_index,
                                  std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_dense_multiply8(std::uint32_t* scratch,std::uint32_t destination,
                                  std::uint32_t a_index,std::uint32_t b_index,
                                  std::uint32_t mask) noexcept {
    if(sirius_native_v4_supported()) sirius_native_v4_dense_multiply8(scratch,destination,a_index,b_index,mask);
    else sirius_native_v3_dense_multiply8(scratch,destination,a_index,b_index,mask);
}
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v3_dense_divide8(std::uint32_t* scratch,std::uint32_t destination,
                                 std::uint32_t a_index,std::uint32_t b_index,
                                 std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v4_dense_divide8(std::uint32_t* scratch,std::uint32_t destination,
                                 std::uint32_t a_index,std::uint32_t b_index,
                                 std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_dense_divide8(std::uint32_t* scratch,std::uint32_t destination,
                                 std::uint32_t a_index,std::uint32_t b_index,
                                 std::uint32_t mask) noexcept {
    if(sirius_native_v4_supported()) sirius_native_v4_dense_divide8(scratch,destination,a_index,b_index,mask);
    else sirius_native_v3_dense_divide8(scratch,destination,a_index,b_index,mask);
}
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v3_dense_sqrt8(std::uint32_t* scratch,std::uint32_t destination,
                              std::uint32_t a_index,std::uint32_t b_index,
                              std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v4_dense_sqrt8(std::uint32_t* scratch,std::uint32_t destination,
                              std::uint32_t a_index,std::uint32_t b_index,
                              std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_dense_sqrt8(std::uint32_t* scratch,std::uint32_t destination,
                              std::uint32_t a_index,std::uint32_t b_index,
                              std::uint32_t mask) noexcept {
    if(sirius_native_v4_supported()) sirius_native_v4_dense_sqrt8(scratch,destination,a_index,b_index,mask);
    else sirius_native_v3_dense_sqrt8(scratch,destination,a_index,b_index,mask);
}
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v3_dense_subtract8(std::uint32_t* scratch,std::uint32_t destination,
                                   std::uint32_t a_index,std::uint32_t b_index,
                                   std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v4_dense_subtract8(std::uint32_t* scratch,std::uint32_t destination,
                                   std::uint32_t a_index,std::uint32_t b_index,
                                   std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_dense_subtract8(std::uint32_t* scratch,std::uint32_t destination,
                                   std::uint32_t a_index,std::uint32_t b_index,
                                   std::uint32_t mask) noexcept {
    if(sirius_native_v4_supported()) sirius_native_v4_dense_subtract8(scratch,destination,a_index,b_index,mask);
    else sirius_native_v3_dense_subtract8(scratch,destination,a_index,b_index,mask);
}
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v3_tagged_add8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v4_tagged_add8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_tagged_add8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept {
    if(sirius_native_v4_supported()) sirius_native_v4_tagged_add8(scratch,tags,destination,a_index,b_index,mask);
    else sirius_native_v3_tagged_add8(scratch,tags,destination,a_index,b_index,mask);
}
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v3_tagged_multiply8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v4_tagged_multiply8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_tagged_multiply8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept {
    if(sirius_native_v4_supported()) sirius_native_v4_tagged_multiply8(scratch,tags,destination,a_index,b_index,mask);
    else sirius_native_v3_tagged_multiply8(scratch,tags,destination,a_index,b_index,mask);
}
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v3_tagged_subtract8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v4_tagged_subtract8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_tagged_subtract8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept {
    if(sirius_native_v4_supported()) sirius_native_v4_tagged_subtract8(scratch,tags,destination,a_index,b_index,mask);
    else sirius_native_v3_tagged_subtract8(scratch,tags,destination,a_index,b_index,mask);
}
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v3_tagged_divide8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v4_tagged_divide8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_tagged_divide8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept {
    if(sirius_native_v4_supported()) sirius_native_v4_tagged_divide8(scratch,tags,destination,a_index,b_index,mask);
    else sirius_native_v3_tagged_divide8(scratch,tags,destination,a_index,b_index,mask);
}
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v3_tagged_sqrt8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_v4_tagged_sqrt8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept;
extern "C" __attribute__((visibility("hidden"))) void sirius_native_tagged_sqrt8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept {
    if(sirius_native_v4_supported()) sirius_native_v4_tagged_sqrt8(scratch,tags,destination,a_index,b_index,mask);
    else sirius_native_v3_tagged_sqrt8(scratch,tags,destination,a_index,b_index,mask);
}
