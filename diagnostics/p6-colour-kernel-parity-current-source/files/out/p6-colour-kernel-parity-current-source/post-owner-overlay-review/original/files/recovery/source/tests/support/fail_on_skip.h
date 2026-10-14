#pragma once

// CMake's GoogleTest discovery treats GTEST_SKIP() as a successful CTest skip.
// The required operational profile force-includes this file in every first-
// party test translation unit so missing capabilities and degenerate fixtures
// become failures. Portable profiles do not include it and retain explicit
// capability-domain skips.

#include <gtest/gtest.h>

#undef GTEST_SKIP
#define GTEST_SKIP() FAIL()
