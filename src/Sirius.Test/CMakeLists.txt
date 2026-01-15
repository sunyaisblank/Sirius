# =============================================================================
# Sirius.Test CMake Configuration
# Google Test-based testing framework
#
# Governance: docs/specification.md
# Test Categories:
#   - Mandatory: Must pass for build to succeed
#   - Correctness: Mathematical invariant tests
#   - Stability: Numerical stability tests
#   - Performance: FPS and resource metrics (soft gate)
# =============================================================================

cmake_minimum_required(VERSION 3.16)
project(SiriusTest LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# === Fetch Google Test ===
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.14.0
)
# Prevent GoogleTest from overriding compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# === Test Executable ===
# NOTE: Legacy metrics (PHMT001A-004A) replaced by unified metric families (PHMT100A-102A)
add_executable(SiriusTests
    # ==========================================================================
    # Unit Tests - Math (TSMT###A)
    # ==========================================================================
    Unit/TSMT001A.cpp    # Dual number arithmetic and derivatives
    Unit/TSMT002A.cpp    # Tensor operations and contractions
    
    # ==========================================================================
    # Unit Tests - Physics (TSPH###A)
    # ==========================================================================
    Unit/TSPH001A.cpp    # Schwarzschild metric properties
    Unit/TSPH002A.cpp    # Christoffel symbol symmetry
    Unit/TSPH003A.cpp    # Geodesic equation integration
    Unit/TSPH004A.cpp    # Kerr metric properties
    Unit/TSPH005A.cpp    # Kerr-Schild coordinates
    Unit/TSPH006A.cpp    # Ellis Drainhole (wormhole)
    Unit/TSPH007A.cpp    # Reissner-Nordström (charged black hole)
    Unit/TSPH008A.cpp    # Numerical metric validation (ADM)
    Unit/TSPH009A.cpp    # RK45 integrator validation
    Unit/TSPH010A.cpp    # Numerical geodesic integrator
    # Unit/TSNR001A.cpp    # [DELETED] NR data loader (obsolete)
    Unit/TSPH011A.cpp    # Spectral rendering utilities
    Unit/TSPH012A.cpp    # Metric derivative validation
    Unit/TSPH014A.cpp    # Double-precision metric tests (Phase 1: Offline Renderer)
    Unit/TSSI001A.cpp    # Symplectic integrator tests (Phase 2: Offline Renderer)
    Unit/TSGD001A.cpp    # Beam propagation tests (Phase 3: Offline Renderer)
    Unit/TSSP001A.cpp    # Spectral radiance tests (Phase 4: Offline Renderer)
    Unit/TSDK001A.cpp    # Accretion disk tests (Phase 5: Offline Renderer)
    # Unit/TSSI002A.cpp    # Scene intersection/anti-aliasing tests (Phase 3-5) [DISABLED: KNAA001A not implemented]
    Unit/TSGD002A.cpp    # Einstein ring magnification tests (Gate G3)
    # Unit/TSSP002A.cpp    # Redshift energy conservation and spectral transport tests [DISABLED: OFTM001A not implemented]
    # === OFFLINE RENDERING TESTS (Production Plan) ===
    # [DISABLED: Sirius.Offline module not yet implemented]
    # Unit/TSOF001A.cpp    # Render session and camera model tests (Phase 6)
    # Unit/TSOF002A.cpp    # Cancellation token tests (Production Plan Phase 1)
    # Unit/TSOF003A.cpp    # Error accumulator tests (Production Plan Phase 1)
    # Unit/TSOF004A.cpp    # Tile scheduler tests (Production Plan Phase 2)
    # Unit/TSOF005A.cpp    # ETA calculator tests (Production Plan Phase 2)
    # Unit/TSOF006A.cpp    # Output driver tests (Production Plan Phase 3)
    # Unit/TSOF007A.cpp    # Checkpoint tests (Production Plan Phase 4)
    # Unit/TSOF008A.cpp    # Display driver tests (Production Plan Phase 5)
    # Unit/TSMB001A.cpp    # Motion blur and temporal sampling tests (Phase 7) [DISABLED: KNST001A not implemented]
    Unit/TSVAL001A.cpp   # DNGR validation and analytic tests (Phase 8)
    
    # ==========================================================================
    # Unit Tests - Camera/DNGR Validation (TSCM###A)
    # ==========================================================================
    Unit/TSCM005A.cpp    # FIDO frame consistency (DNGR Eq A.2-A.7)
    Unit/TSCM006A.cpp    # G-factor validation (DNGR Eq A.16)
    Unit/TSCM007A.cpp    # Geodesic deviation (DNGR Eq A.18-A.27)
    
    # ==========================================================================
    # Diagnostic Tests (TSDG###A) - Numerical Stability
    # ==========================================================================

    Diagnostic/TSDG001A.cpp  # Edge cases and precision
    Diagnostic/TSDG002A.cpp  # NaN/Inf detection
    Diagnostic/TSDG003A.cpp  # Conservation law verification
    Diagnostic/TSDG004A.cpp  # Determinism verification
    
    # ==========================================================================
    # Integration Tests (TSIN###A) - End-to-end validation
    # ==========================================================================
    Integration/TSIN001A.cpp  # Geodesic path validation
    Integration/TSIN002A.cpp  # Metric loader chain (Unified Families)
    Integration/TSIN003A.cpp  # Render pipeline integration
    # Integration/TSIN004A.cpp  # Offline renderer integration (Production Plan) [DISABLED: Sirius.Offline not implemented]
    
    # ==========================================================================
    # Benchmark Tests (TSBM###A) - Performance validation
    # ==========================================================================
    Benchmark/TSBM001A.cpp   # Integration accuracy (Killing vectors)
    Benchmark/TSBM002A.cpp   # FPS performance tracking
    Benchmark/TSBM003A.cpp   # Christoffel computation: Spherical vs Cartesian
    # Benchmark/TSBM004A.cpp   # [ARCHIVED] Symplectic vs RK45 integrator comparison
    
    # ==========================================================================
    # Performance Tests (TSPF###A) - Resource usage validation
    # ==========================================================================
    Performance/TSPF001A.cpp  # FPS thresholds
    Performance/TSPF002A.cpp  # Memory usage
)

# === Include Directories ===
target_include_directories(SiriusTests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Tensor
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Autodiff
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Metric
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Geodesic
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Coordinate
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Spectral
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Disk
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Symplectic
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Transport
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Render/Session
    ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/glm
)

# === Link Libraries ===
target_link_libraries(SiriusTests PRIVATE
    SiriusCore
    GTest::gtest_main
)

# === Enable CTest ===
include(GoogleTest)
gtest_discover_tests(SiriusTests
    PROPERTIES
        LABELS "Sirius"
)

# === Output Directory ===
set_target_properties(SiriusTests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
)

# =============================================================================
# TEST LABELING FOR MANDATORY GATE
# =============================================================================
# After tests are discovered, we apply labels via test properties.
# GoogleTest discovery happens at build time, so we use PROPERTIES in discover.
#
# Label Definitions:
#   Mandatory   - Build fails if these tests fail
#   Correctness - Mathematical invariant tests  
#   Stability   - Numerical stability tests
#   Performance - FPS/resource metrics (warning only)
# =============================================================================

# Note: gtest_discover_tests applies labels at configure time.
# For fine-grained control, individual test labels can be set via:
#   set_tests_properties(TestName PROPERTIES LABELS "Mandatory;Correctness")
#
# Since GoogleTest creates tests dynamically, we use a post-discovery approach:
# The following regex patterns match test names to apply labels.

# After build, label tests using CTest regex matching
# This is executed via a custom script or manually after first configure

message(STATUS "")
message(STATUS "=== Sirius Test Suite ===")
message(STATUS "  Unit Tests:        13 files (TSMT, TSPH, TSNR)")
message(STATUS "  Diagnostic Tests:  4 files (TSDG)")
message(STATUS "  Integration Tests: 3 files (TSIN)")
message(STATUS "  Benchmark Tests:   2 files (TSBM)")
message(STATUS "  Performance Tests: 2 files (TSPF)")
message(STATUS "")
message(STATUS "  Run all tests:     ctest --test-dir <build>")
message(STATUS "  Run mandatory:     ctest --test-dir <build> -L Mandatory")
message(STATUS "")

