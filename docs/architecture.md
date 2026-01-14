# Sirius Architecture

*Structural Design and Component Organisation*

---

## Overview

This document describes the structural organisation of Sirius, including architectural principles, layer relationships, component naming conventions, and the component registry. It provides the comprehensive reference for understanding where code lives and how it is identified.

---

## Part I: Architectural Principles

### 1.1 Separation of Concerns

Each component has a single, well-defined responsibility:

| Layer | Responsibility |
|-------|----------------|
| Physics | Metric tensors, geodesics, spacetime geometry |
| Mathematics | Tensors, dual numbers, linear algebra |
| Rendering | Ray tracing, shaders, display |
| Application | Window management, UI, orchestration |
| Testing | Validation, benchmarks, diagnostics |

### 1.2 Dependency Direction

Dependencies flow in a controlled direction:

```
Application → Rendering → Physics → Mathematics
          ↘    UI    ↗
```

Lower layers have no knowledge of higher layers. Physics does not know about rendering; mathematics does not know about physics.

### 1.3 Interface Boundaries

Components communicate through explicit interfaces. The `IMetric` interface allows the renderer to work with any spacetime geometry without coupling to specific implementations.

---

## Part II: Layer Diagram

```mermaid
flowchart TB
    subgraph INFRASTRUCTURE["INFRASTRUCTURE LAYER"]
        Infra["<b>Sirius.Infrastructure</b><br/>Window, UI, plugins,<br/>application orchestration"]
    end

    subgraph RENDER["RENDERING LAYER"]
        OptiX["<b>Sirius.Render/OptiX</b><br/>RTX-accelerated raytracing"]
        Session["<b>Sirius.Render/Session</b><br/>Batch rendering,<br/>session management"]
        Pipeline["<b>Sirius.Render/Pipeline</b><br/>Display, tonemapping"]
    end

    subgraph CORE["CORE LAYER (Domain Primitives)"]
        Metric["<b>Sirius.Core/Metric</b><br/>Metric tensors"]
        Geodesic["<b>Sirius.Core/Geodesic</b><br/>Geodesic integration"]
        Camera["<b>Sirius.Core/Camera</b><br/>Observer motion"]
        Spectral["<b>Sirius.Core/Spectral</b><br/>Radiative transfer"]
        Tensor["<b>Sirius.Core/Tensor</b><br/>Tensor algebra"]
    end

    subgraph TEST["TEST LAYER"]
        Test["<b>Sirius.Test</b><br/>Unit, diagnostic,<br/>benchmark tests"]
    end

    INFRASTRUCTURE --> RENDER
    RENDER --> CORE
    TEST -.->|validates| CORE
    TEST -.->|validates| RENDER

    style INFRASTRUCTURE fill:#e1f5fe
    style RENDER fill:#fff3e0
    style CORE fill:#e8f5e9
    style TEST fill:#fce4ec
```

---

## Part III: Directory Structure

```
Sirius/
├── docs/                         # Documentation
│   ├── README.md                 # Documentation index
│   ├── philosophy.md             # Design philosophy
│   ├── foundations.md            # Mathematical theory
│   ├── types.md                  # Type system
│   ├── specification.md          # Formal requirements
│   ├── guide.md                  # Practical operation
│   ├── standard.md               # Coding conventions
│   ├── architecture.md           # This document
│   └── refactor.md               # Refactoring discipline
│
├── src/                          # Source code
│   ├── Sirius.Core/              # Domain primitives
│   │   ├── Tensor/               # Multi-dimensional tensors
│   │   ├── Autodiff/             # Automatic differentiation
│   │   ├── Metric/               # Metric tensor implementations
│   │   ├── Geodesic/             # Geodesic integration
│   │   ├── Coordinate/           # Coordinate transformations
│   │   ├── Spectral/             # Radiative transfer
│   │   ├── Disk/                 # Accretion disk physics
│   │   ├── Symplectic/           # Symplectic integrators
│   │   ├── Transport/            # Parallel transport
│   │   └── Camera/               # Geodesic camera motion
│   │
│   ├── Sirius.Render/            # Rendering layer
│   │   ├── OptiX/                # RTX-accelerated raytracing
│   │   ├── Shader/               # GLSL shaders
│   │   ├── Pipeline/             # Display pipeline
│   │   ├── Texture/              # Starfield, resources
│   │   ├── Session/              # Batch rendering (offline)
│   │   └── Kernel/               # Compute kernels
│   │
│   ├── Sirius.Infrastructure/    # Application orchestration
│   │   ├── CRAP001A.h/cpp        # Application class
│   │   ├── CRWN001A.h/cpp        # Window management
│   │   ├── CRPM001A.h/cpp        # Plugin manager
│   │   ├── CREP001A.cpp          # Entry point (main)
│   │   └── UIMN001A.h/cpp        # UI manager
│   │
│   └── Sirius.Test/              # Test suite
│       ├── Unit/                 # Unit tests
│       ├── Diagnostic/           # Numerical stability tests
│       └── Benchmark/            # Performance benchmarks
│
├── lib/                          # External dependencies
│   ├── glad/                     # OpenGL loader
│   ├── glfw/                     # Window management
│   ├── glm/                      # Math library
│   ├── imgui/                    # UI framework
│   ├── stb/                      # Image I/O
│   └── tinyexr/                  # EXR format support
│
├── bin/                          # Build output (gitignored)
│   └── Sirius.Build/             # CMake build directory
│
├── ren/                          # Render output (gitignored)
│
├── CMakeLists.txt                # Build configuration
└── README.md                     # Project overview
```

---

## Part IV: Component Naming Convention

### 4.1 Format

```
[Domain][Category][Sequence][Variant]
```

**Example:** `PHMT001A`

- **Domain:** `PH` (Physics)
- **Category:** `MT` (Metric Tensor)
- **Sequence:** `001` (First component)
- **Variant:** `A` (Primary implementation)

### 4.2 Domain Codes

| Code | Domain | Description |
|------|--------|-------------|
| `PH` | Physics | Metrics, geodesics, spacetime geometry |
| `MT` | Math | Tensors, dual numbers, algebra |
| `RD` | Render | Raytracer, shaders, display |
| `UI` | Interface | User controls, ImGui |
| `CR` | Core | Application, windowing, plugins (Infrastructure) |
| `CM` | Camera | Camera models and motion |
| `OF` | Offline | Batch rendering, sessions |
| `PP` | PostProcess | Tonemapping, bloom, effects |
| `NR` | Numerics | Data loaders, converters |
| `TS` | Test | Test cases |

### 4.3 Category Codes

#### Physics (PH)

| Code | Category | Description |
|------|----------|-------------|
| `MT` | Metric | Metric tensor implementations |
| `GD` | Geodesic | Geodesic integration |
| `CF` | Christoffel | Connection coefficients |
| `CT` | Coordinate | Coordinate transformations |
| `OB` | Observer | Observer physics, tetrads |
| `RS` | Redshift | Frequency shift calculations |

#### Mathematics (MT)

| Code | Category | Description |
|------|----------|-------------|
| `DL` | Dual | Dual numbers, autodiff |
| `TN` | Tensor | Tensor algebra |
| `VT` | Vector | Vector types |
| `QT` | Quaternion | Rotation mathematics |

#### Rendering (RD)

| Code | Category | Description |
|------|----------|-------------|
| `RT` | Raytracer | Core raytracing logic |
| `OP` | OptiX | OptiX host/device code |
| `TX` | Texture | Texture handling |
| `DP` | Display | OpenGL display/window |

#### Test Suite (TS)

| Code | Category | Description |
|------|----------|-------------|
| `UN` | Unit | Component-level unit tests |
| `IN` | Integration | End-to-end integration tests |
| `DG` | Diagnostic | Numerical stability tests |
| `BM` | Benchmark | Performance benchmarks |

### 4.4 Variant Codes

| Code | Meaning | Usage |
|------|---------|-------|
| `A` | Primary | Production-ready implementation |
| `B` | Alternative | Alternative algorithm |
| `X` | Experimental | Under development |
| `D` | Deprecated | Superseded, kept for reference |

---

## Part V: Component Registry

### 5.1 Sirius.Math

| Code | File | Description |
|------|------|-------------|
| `MTDL001A` | `MTDL001A.h` | Dual number arithmetic |
| `MTTN001A` | `MTTN001A.h/cpp` | Multi-dimensional tensor |

### 5.2 Sirius.Physics

| Code | File | Description |
|------|------|-------------|
| `PHMT000A` | `PHMT000A.h` | Metric tensor interface |
| `PHMT001A` | `PHMT001A.h/cpp` | Minkowski metric |
| `PHMT002A` | `PHMT002A.h/cpp` | Schwarzschild metric |
| `PHMT003A` | `PHMT003A.h/cpp` | Kerr metric |
| `PHMT004A` | `PHMT004A.h/cpp` | Reissner-Nordström metric |
| `PHMT100A` | `PHMT100A.h` | Kerr-Schild family (unified) |
| `PHMT101A` | `PHMT101A.h` | Morris-Thorne family |
| `PHMT102A` | `PHMT102A.h` | Warp drive family |
| `PHGD001A` | `PHGD001A.h/cpp` | Geodesic integrator |
| `PHGD002A` | `PHGD002A.h/cpp` | Numerical geodesic integrator |
| `PHCT001A` | `PHCT001A.h` | Coordinate transformations |

### 5.3 Sirius.Render

| Code | File | Description |
|------|------|-------------|
| `RDRT001A` | `RDRT001A.cpp` | Main renderer |
| `RDOP001A` | `RDOP001A.cu` | OptiX host code |
| `RDOP002A` | `RDOP002A.cu` | OptiX device kernels |
| `RDOP003A` | `RDOP003A.h` | Launch parameters |

### 5.4 Sirius.Core

| Code | File | Description |
|------|------|-------------|
| `CREP001A` | `CREP001A.cpp` | Entry point (main) |
| `CRAP001A` | `CRAP001A.h/cpp` | Application class |
| `CRPM001A` | `CRPM001A.h/cpp` | Plugin manager |
| `CRWN001A` | `CRWN001A.h/cpp` | Window management |

### 5.5 Sirius.Offline

| Code | File | Description |
|------|------|-------------|
| `OFRS001A` | `OFRS001A.h` | Render session |
| `OFCM001A` | `OFCM001A.h` | Camera model |
| `OFIO001A` | `OFIO001A.h` | Image I/O |
| `OFCT001A` | `OFCT001A.h` | Cancellation token |
| `OFER001A` | `OFER001A.h` | Error accumulator |

### 5.6 Sirius.Test

| Code | File | Description |
|------|------|-------------|
| `TSMT001A` | `Unit/TSMT001A.cpp` | Dual number tests |
| `TSMT002A` | `Unit/TSMT002A.cpp` | Tensor operation tests |
| `TSPH001A` | `Unit/TSPH001A.cpp` | Schwarzschild metric tests |
| `TSPH002A` | `Unit/TSPH002A.cpp` | Christoffel symbol tests |
| `TSPH003A` | `Unit/TSPH003A.cpp` | Geodesic equation tests |
| `TSDG001A` | `Diagnostic/TSDG001A.cpp` | Numerical stability tests |
| `TSBM001A` | `Benchmark/TSBM001A.cpp` | Integration accuracy tests |

---

## Part VI: Data Flow

### 6.1 Ray Tracing Pipeline

```mermaid
flowchart LR
    Camera["Camera"] --> InitRay["Initialise Ray"]
    InitRay --> Integrate["Geodesic Integration"]
    Integrate --> Terminate{"Terminated?"}
    Terminate -->|No| Metric["Evaluate Metric"]
    Metric --> Christoffel["Compute Christoffel"]
    Christoffel --> Step["RK4/DOPRI Step"]
    Step --> Integrate
    Terminate -->|Escaped| Background["Sample Background"]
    Terminate -->|Disk| Emission["Compute Emission"]
    Terminate -->|Horizon| Black["Return Black"]
    Background --> Pixel["Write Pixel"]
    Emission --> Pixel
    Black --> Pixel
```

### 6.2 Metric Evaluation Flow

```mermaid
flowchart TB
    Position["Position (x,y,z,t)"] --> Validate["Validate Domain"]
    Validate -->|Valid| Compute["Compute Metric Components"]
    Validate -->|Invalid| Clamp["Clamp to Boundary"]
    Clamp --> Compute
    Compute --> Inverse["Compute Inverse"]
    Inverse --> Christoffel["Compute Christoffel"]
    Christoffel --> Return["Return State"]
```

---

## Part VII: File Naming

### Source Files

- Headers: `ComponentCode.h`
- Implementation: `ComponentCode.cpp`
- CUDA: `ComponentCode.cu`

### Shaders

- Compute: `*.comp`
- GLSL includes: `*.glsl`
- Vertex: `*.vert`
- Fragment: `*.frag`

---

## Appendix A: Compliance

All changes to the codebase MUST adhere to this standard. Non-compliant components will be flagged during code review. New components MUST be registered in the component registry before merging.

---

*End of Architecture*
