// Separately named strict copy; original GCC/public math remains untouched.
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <immintrin.h>
#include <limits>
#include <initializer_list>
#pragma STDC FENV_ACCESS ON
#pragma clang fp exceptions(strict)
#pragma clang fp contract(off)
#define sirius_native_dense_add8 sirius_scope_v4_dense_add8
#define sirius_native_dense_multiply8 sirius_scope_v4_dense_multiply8
#define sirius_native_dense_divide8 sirius_scope_v4_dense_divide8
#define sirius_native_dense_sqrt8 sirius_scope_v4_dense_sqrt8
#define sirius_native_dense_subtract8 sirius_scope_v4_dense_subtract8
#define sirius_native_tagged_add8 sirius_scope_v4_tagged_add8
#define sirius_native_tagged_multiply8 sirius_scope_v4_tagged_multiply8
#define sirius_native_tagged_subtract8 sirius_scope_v4_tagged_subtract8
#define sirius_native_tagged_divide8 sirius_scope_v4_tagged_divide8
#define sirius_native_tagged_sqrt8 sirius_scope_v4_tagged_sqrt8
#include "/home/astra/.project/Sirius/out/native-cpu-exact-product-zero-residual/retained8.cpp"
