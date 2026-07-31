#pragma once

// Compile-time language capability facts. These names distinguish native
// standard support from Sirius's enforced substitutes; neither the build nor
// runtime reporting may infer P2900/P2996 support from C++26 mode alone.

namespace sirius::base {

#if defined(__cpp_contracts) && __cpp_contracts >= 202502L
inline constexpr bool kNativeP2900Contracts = true;
#else
inline constexpr bool kNativeP2900Contracts = false;
#endif

#if (defined(__cpp_impl_reflection) && __cpp_impl_reflection >= 202506L) || \
    (defined(__cpp_reflection) && __cpp_reflection >= 202506L)
inline constexpr bool kNativeP2996Reflection = true;
#else
inline constexpr bool kNativeP2996Reflection = false;
#endif

enum class ContractImplementation {
    NativeP2900,
    CheckedMacro,
};

enum class ReflectionImplementation {
    NativeP2996,
    ExplicitSchema,
};

inline constexpr ContractImplementation kContractImplementation =
    kNativeP2900Contracts ? ContractImplementation::NativeP2900
                          : ContractImplementation::CheckedMacro;

// Sirius currently binds configuration explicitly with nlohmann serializers
// and packs kernel parameters explicitly. When P2996 is available, adopting it
// requires a deliberate implementation change; a feature macro alone must not
// cause an unreviewed representation change.
inline constexpr ReflectionImplementation kReflectionImplementation =
    ReflectionImplementation::ExplicitSchema;

}  // namespace sirius::base
