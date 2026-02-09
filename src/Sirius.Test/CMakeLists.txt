# =============================================================================
# Sirius.Test CMake Configuration
# Google Test-based testing framework
#
# Governance: docs/specification.md
# Test Categories (via CTest labels):
#   - Mandatory: Must pass for build to succeed
#   - Correctness: Mathematical invariant tests
#   - Stability: Numerical stability tests
#   - Performance: FPS and resource metrics (soft gate)
#
# Directory structure mirrors src/ for 1:1 coverage mapping.
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
add_executable(SiriusTests
    # ==========================================================================
    # Sirius.Core — Autodiff
    # ==========================================================================
    Sirius.Core/Autodiff/TSMT001A.cpp    # Dual number arithmetic and derivatives

    # ==========================================================================
    # Sirius.Core — Tensor
    # ==========================================================================
    Sirius.Core/Tensor/TSMT002A.cpp      # Tensor operations and contractions

    # ==========================================================================
    # Sirius.Core — Metric
    # ==========================================================================
    Sirius.Core/Metric/TSPH001A.cpp      # Schwarzschild metric properties
    Sirius.Core/Metric/TSPH002A.cpp      # Christoffel symbol symmetry
    Sirius.Core/Metric/TSPH004A.cpp      # Kerr metric properties
    Sirius.Core/Metric/TSPH005A.cpp      # Kerr-Schild coordinates
    Sirius.Core/Metric/TSPH006A.cpp      # Ellis Drainhole (wormhole)
    Sirius.Core/Metric/TSPH007A.cpp      # Reissner-Nordstrom (charged black hole)
    Sirius.Core/Metric/TSPH008A.cpp      # Numerical metric validation (ADM)
    Sirius.Core/Metric/TSPH012A.cpp      # Metric derivative validation
    Sirius.Core/Metric/TSPH014A.cpp      # Double-precision metric tests
    Sirius.Core/Metric/TSCM005A.cpp      # FIDO frame consistency (DNGR Eq A.2-A.7)
    Sirius.Core/Metric/TSCM006A.cpp      # G-factor validation (DNGR Eq A.16)
    Sirius.Core/Metric/TSCM007A.cpp      # Geodesic deviation (DNGR Eq A.18-A.27)
    Sirius.Core/Metric/TSSM001A.cpp      # SMBH parameter scaling and invariants
    Sirius.Core/Metric/TSVAL001A.cpp     # DNGR validation and analytic tests
    Sirius.Core/Metric/TSDG007A.cpp      # Comprehensive metric tensor validation
    Sirius.Core/Metric/TSDG008A.cpp      # Analytic solution validation (photon sphere, ISCO, deflection)
    Sirius.Core/Metric/TSIN002A.cpp      # Metric loader chain (Unified Families)
    Sirius.Core/Metric/TSMT102A.cpp      # Alcubierre warp drive metric properties [NEW]

    # ==========================================================================
    # Sirius.Core — Geodesic
    # ==========================================================================
    Sirius.Core/Geodesic/TSPH003A.cpp    # Geodesic equation integration
    Sirius.Core/Geodesic/TSPH009A.cpp    # RK45 integrator validation
    Sirius.Core/Geodesic/TSPH010A.cpp    # Numerical geodesic integrator
    Sirius.Core/Geodesic/TSGD001A.cpp    # Beam propagation tests
    Sirius.Core/Geodesic/TSGD002A.cpp    # Einstein ring magnification tests
    Sirius.Core/Geodesic/TSDG001A.cpp    # Edge cases and precision
    Sirius.Core/Geodesic/TSDG002A.cpp    # NaN/Inf detection
    Sirius.Core/Geodesic/TSDG003A.cpp    # Conservation law verification
    Sirius.Core/Geodesic/TSDG004A.cpp    # Determinism verification
    Sirius.Core/Geodesic/TSDG006A.cpp    # Mandatory Killing vector tests
    Sirius.Core/Geodesic/TSBM001A.cpp    # Integration accuracy (Killing vectors)
    Sirius.Core/Geodesic/TSIN001A.cpp    # Geodesic path validation

    # ==========================================================================
    # Sirius.Core — Disk
    # ==========================================================================
    Sirius.Core/Disk/TSDK001A.cpp        # Accretion disk tests
    Sirius.Core/Disk/TSDK002A.cpp        # Volumetric disk tests
    Sirius.Core/Disk/TSTR001A.cpp        # Turbulence model power spectrum validation
    Sirius.Core/Disk/TSCR001A.cpp        # Corona model tests [NEW]

    # ==========================================================================
    # Sirius.Core — Jet
    # ==========================================================================
    Sirius.Core/Jet/TSJT002A.cpp         # MHD jet magnetic flux conservation
    Sirius.Core/Jet/TSJT001A.cpp         # Relativistic jet emission model [NEW]

    # ==========================================================================
    # Sirius.Core — Spectral
    # ==========================================================================
    Sirius.Core/Spectral/TSPH011A.cpp    # Spectral rendering utilities
    Sirius.Core/Spectral/TSSP001A.cpp    # Spectral radiance tests
    Sirius.Core/Spectral/TSSP002A.cpp    # Spectral reference validation (Planck, Wien, Stefan-Boltzmann)

    # ==========================================================================
    # Sirius.Core — Coordinate
    # ==========================================================================
    Sirius.Core/Coordinate/TSPH015A.cpp  # Coordinate transformation tests (PHCT002A round-trip)
    Sirius.Core/Coordinate/TSCT001A.cpp  # Coordinate transform interface tests [NEW]

    # ==========================================================================
    # Sirius.Core — Camera
    # ==========================================================================
    Sirius.Core/Camera/TSCM001A.cpp      # Camera lens model tests

    # ==========================================================================
    # Sirius.Core — Symplectic
    # ==========================================================================
    Sirius.Core/Symplectic/TSSI001A.cpp  # Symplectic integrator tests
    Sirius.Core/Symplectic/TSSI002A.cpp  # Beam sampling and scene intersection tests

    # ==========================================================================
    # Sirius.Core — Polarisation
    # ==========================================================================
    Sirius.Core/Polarisation/TSPL001A.cpp  # Stokes vector and Mueller matrix tests [NEW]

    # ==========================================================================
    # Sirius.Core — Transport
    # ==========================================================================
    Sirius.Core/Transport/TSTP001A.cpp   # Vec4d, GeodesicStateD, Mat4d tests [NEW]

    # ==========================================================================
    # Sirius.Core — PostProcess
    # ==========================================================================
    Sirius.Core/PostProcess/TSPP001A.cpp # Tonemapping operator tests [NEW]

    # ==========================================================================
    # Sirius.Core — Environment
    # ==========================================================================
    Sirius.Core/Environment/TSSF001A.cpp # Starfield generator tests [NEW]

    # ==========================================================================
    # Sirius.Render — Session
    # ==========================================================================
    Sirius.Render/Session/TSOF001A.cpp   # Render session tests
    Sirius.Render/Session/TSOF002A.cpp   # Cancellation token tests
    Sirius.Render/Session/TSOF003A.cpp   # Error accumulator tests
    Sirius.Render/Session/TSOF004A.cpp   # Tile scheduler tests
    Sirius.Render/Session/TSOF005A.cpp   # ETA calculator tests
    Sirius.Render/Session/TSOF006A.cpp   # Output driver tests
    Sirius.Render/Session/TSOF007A.cpp   # Checkpoint tests
    Sirius.Render/Session/TSOF008A.cpp   # Display driver tests
    Sirius.Render/Session/TSOF009A.cpp   # PNG output driver tests
    Sirius.Render/Session/TSOF010A.cpp   # Offline renderer integration tests
    Sirius.Render/Session/TSMB001A.cpp   # Motion blur and temporal sampling tests

    # ==========================================================================
    # Sirius.Render — Integration
    # ==========================================================================
    Sirius.Render/Integration/TSIN005A.cpp  # GeodesicTracer integration tests
    Sirius.Render/Integration/TSDG005A.cpp  # GPU conservation validation
    Sirius.Render/Integration/TSIN006A.cpp  # GPU/CPU parity tests

    # ==========================================================================
    # Sirius.Render — Pipeline
    # ==========================================================================
    Sirius.Render/Pipeline/TSIN003A.cpp  # Render pipeline integration

    # ==========================================================================
    # Sirius.Render — Offline
    # ==========================================================================
    Sirius.Render/Offline/TSIN004A.cpp   # Offline renderer integration

    # ==========================================================================
    # Sirius.Infrastructure — Configuration
    # ==========================================================================
    Sirius.Infrastructure/Configuration/TSFL001A.cpp  # Film simulation unit tests

    # ==========================================================================
    # Sirius.Core — Metric (Benchmark & Performance)
    # ==========================================================================
    Sirius.Core/Metric/TSBM002A.cpp      # FPS performance tracking
    Sirius.Core/Metric/TSBM003A.cpp      # Christoffel computation: Spherical vs Cartesian
    Sirius.Core/Metric/TSPF001A.cpp      # FPS thresholds
    Sirius.Core/Metric/TSPF002A.cpp      # Memory usage

    # ==========================================================================
    # Implementation sources (required by tests)
    # ==========================================================================
    ../Sirius.Render/Output/OUPN001A.cpp   # PNG writer implementation
    ../Sirius.Render/Output/OUEW001A.cpp   # EXR writer implementation
    ../../lib/tinyexr/miniz.c              # miniz compression (required by tinyexr)
    ../Sirius.Render/Integration/GTRC001A.cpp  # GeodesicTracer implementation
)

# === Include Directories ===
# Note: SiriusCore PUBLIC includes are inherited via target_link_libraries
target_include_directories(SiriusTests PRIVATE
    # Test fixtures and infrastructure (TSFT*.h at test root)
    ${CMAKE_CURRENT_SOURCE_DIR}
    # Source root (for "Sirius.Render/..." style includes)
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    # Render module directories
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Render
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Render/Integration
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Render/Acceleration
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Render/Acceleration/OptiX
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Render/Session
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Render/Scheduling
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Render/Output
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Render/Camera
    # Infrastructure module directories
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Infrastructure
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Infrastructure/Configuration
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Infrastructure/Session
    # Core module directories
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Disk
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Jet
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Metric
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Environment
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Polarisation
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Transport
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/PostProcess
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Camera
    ${CMAKE_CURRENT_SOURCE_DIR}/../Sirius.Core/Coordinate
    # Third-party image I/O libraries
    ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/stb
    ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/tinyexr
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
# Labels are applied by test fixture name in CTestLabels.cmake.
# Directory structure encodes module affinity; labels encode test category.
# Included via TEST_INCLUDE_FILES so labels apply after gtest discovery.
# =============================================================================
set_property(DIRECTORY APPEND PROPERTY
    TEST_INCLUDE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/CTestLabels.cmake"
)

message(STATUS "")
message(STATUS "=== Sirius Test Suite ===")
message(STATUS "  Sirius.Core:           59 files (Autodiff, Tensor, Metric, Geodesic, Disk, Jet, Spectral, Coordinate, Symplectic, Camera, Polarisation, Transport, PostProcess, Environment)")
message(STATUS "  Sirius.Render:         15 files (Session, Integration, Pipeline, Offline)")
message(STATUS "  Sirius.Infrastructure:  1 file  (Configuration)")
message(STATUS "")
message(STATUS "  Run all tests:     ctest --test-dir <build>")
message(STATUS "  Run mandatory:     ctest --test-dir <build> -L Mandatory")
message(STATUS "")
