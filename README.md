# Sirius

Sirius is a general relativistic ray tracing engine. It traces photon geodesics through analytically-specified metric tensor fields, rendering gravitational lensing, black hole shadows, and accretion disk emission near compact objects.

## Geodesic Ray Tracing

Sirius solves the geodesic equation of general relativity, the equation governing the path of a photon through curved spacetime:

$$
\frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\nu\rho} \frac{dx^\nu}{d\lambda} \frac{dx^\rho}{d\lambda} = 0
$$

where $x^\mu$ denotes the spacetime coordinates, $\lambda$ is an affine parameter along the ray, and $\Gamma^\mu_{\nu\rho}$ are the Christoffel symbols encoding spacetime curvature. The Christoffel symbols are computed analytically from the metric tensor:

$$
\Gamma^\mu_{\nu\rho} = \frac{1}{2} g^{\mu\sigma} \left( \partial_\nu g_{\sigma\rho} + \partial_\rho g_{\sigma\nu} - \partial_\sigma g_{\nu\rho} \right)
$$

Sirius computes these connection coefficients analytically via dual number automatic differentiation, eliminating the truncation error that numerical differentiation would introduce across hundreds of integration steps per ray. See [Design Rationale](docs/philosophy.md) for why this matters.

## Metric Tensor Families

Spacetime geometries are organised into parameterised families. The Kerr-Schild family encompasses all stationary black hole solutions through a four-parameter space $(M, a, Q, \Lambda)$: mass, angular momentum, electric charge, and cosmological constant. The metric takes the Kerr-Schild form:

$$
g_{\mu\nu} = \eta_{\mu\nu} + H \, l_\mu l_\nu
$$

where $\eta_{\mu\nu}$ is the Minkowski metric, $H$ is a scalar function encoding gravitational field strength, and $l_\mu$ is a null vector field. This Cartesian representation yields polynomial Christoffel symbols, avoiding the coordinate singularities at the rotation poles that plague Boyer-Lindquist coordinates.

The Morris-Thorne family covers traversable wormhole geometries. The Alcubierre family implements warp drive metrics.

## Radiative Transfer

Beyond geodesic tracing, the system implements relativistic radiative transfer for volumetric emission and absorption along the ray path:

$$
\frac{dI_\nu}{d\lambda} = j_\nu - \alpha_\nu I_\nu
$$

where $I_\nu$ is the specific intensity, $j_\nu$ is the emission coefficient, and $\alpha_\nu$ is the absorption coefficient. This enables rendering of accretion disks with physically-motivated Novikov-Thorne temperature profiles and Doppler beaming.

## Gravitational Lensing

Strong lensing near black holes creates multiple images and Einstein rings when a background source aligns with the gravitational lens. The photon sphere at $r = 3M$ (for Schwarzschild geometry) is the surface of unstable circular photon orbits, producing the characteristic black hole shadow.

## Architecture

```
Application (Sirius.Infrastructure) → Rendering (Sirius.Render) → Physics (Sirius.Core)
```

Dependencies flow downward. Lower layers have no knowledge of higher layers.

| Component | Responsibility |
|-----------|----------------|
| `Sirius.Core` | Domain layer: tensors, metrics, geodesics, coordinate transforms |
| `Sirius.Render` | GPU acceleration: OptiX ray tracing, shaders, render pipeline |
| `Sirius.Infrastructure` | Application services: windowing, UI, configuration |
| `Sirius.Test` | Unit, diagnostic, integration, and benchmark tests |

## Project Structure

```
Sirius/
├── docs/                    # Documentation
├── lib/                     # External dependencies
│   ├── glad/                # OpenGL loader
│   ├── glfw/                # Windowing
│   ├── glm/                 # Mathematics
│   ├── imgui/               # Immediate mode GUI
│   ├── stb/                 # Image loading
│   └── tinyexr/             # EXR output
├── src/                     # Source code
│   ├── Sirius.Core/         # Domain layer
│   │   ├── Tensor/          # Tensor algebra
│   │   ├── Autodiff/        # Dual number differentiation
│   │   ├── Metric/          # Metric tensor families
│   │   ├── Geodesic/        # Geodesic integration
│   │   ├── Coordinate/      # Coordinate transforms
│   │   ├── Spectral/        # Spectral rendering
│   │   ├── Disk/            # Accretion disk physics
│   │   ├── Symplectic/      # Symplectic integrators
│   │   └── Transport/       # Radiative transfer
│   ├── Sirius.Render/       # Rendering layer
│   │   ├── Acceleration/    # OptiX and GPU backends
│   │   ├── Session/         # Batch rendering
│   │   ├── Pipeline/        # Render pipeline
│   │   ├── Camera/          # Camera models
│   │   ├── Output/          # EXR, PNG output
│   │   └── Buffer/          # Frame buffers
│   ├── Sirius.Infrastructure/ # Application services
│   └── Sirius.Test/         # Test suite
├── bin/                     # Build output
├── renders/                 # Render output (gitignored)
└── CMakeLists.txt           # Build configuration
```

## Documentation

The `docs/` directory is designed to be read in order:

| # | Document | Covers |
|---|----------|--------|
| 1 | **[Design Rationale](docs/philosophy.md)** | Why the system is built the way it is |
| 2 | **[Foundations](docs/foundations.md)** | Mathematical theory, from differential geometry to numerical integration |
| 3 | **[Types](docs/type.md)** | Type system: primitives, coordinates, tensors, metric interface |
| 4 | **[Specification](docs/specification.md)** | Performance targets, invariant tolerances, precision requirements |
| 5 | **[Guide](docs/guide.md)** | Building, configuring, and running the renderer |
| 6 | **[Standard](docs/standard.md)** | Coding conventions (JPL/MISRA-based) |
| 7 | **[Architecture](docs/architecture.md)** | Layer structure, naming convention, component registry |
| 8 | **[Refactoring](docs/refactor.md)** | Method for changing code safely |
| 9 | **[Commentary](docs/comment.md)** | Inline documentation and comment standards |
