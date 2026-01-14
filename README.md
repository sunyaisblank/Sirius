# Sirius

Sirius is a general relativistic ray tracing engine designed to render visualisations of curved spacetime. The system implements exact geodesic integration on analytically-specified metric tensor fields, enabling the simulation of light propagation near black holes, wormholes, and other exotic geometries predicted by general relativity.

## Geodesic Ray Tracing

The theoretical foundation rests on the geodesic equation of general relativity, which describes the path of massless particles through curved spacetime. A photon's trajectory satisfies the geodesic equation:

$$
\frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\nu\rho} \frac{dx^\nu}{d\lambda} \frac{dx^\rho}{d\lambda} = 0
$$

where $x^\mu$ denotes the spacetime coordinates, $\lambda$ is an affine parameter along the ray, and $\Gamma^\mu_{\nu\rho}$ are the Christoffel symbols encoding the curvature of spacetime. The Christoffel symbols are computed analytically from the metric tensor via:

$$
\Gamma^\mu_{\nu\rho} = \frac{1}{2} g^{\mu\sigma} \left( \partial_\nu g_{\sigma\rho} + \partial_\rho g_{\sigma\nu} - \partial_\sigma g_{\nu\rho} \right)
$$

The system enforces exclusively analytic computation of these connection coefficients, eliminating numerical differentiation errors and enabling exact geodesic integration.

## Metric Tensor Families

The architecture organises spacetime geometries into parameterised families, where continuous variation of physical parameters interpolates between distinct spacetimes. The Kerr-Schild family encompasses all stationary black hole solutions through a four-parameter space $(M, a, Q, \Lambda)$ representing mass, angular momentum, electric charge, and cosmological constant respectively. The metric tensor takes the Kerr-Schild form:

$$
g_{\mu\nu} = \eta_{\mu\nu} + H \, l_\mu l_\nu
$$

where $\eta_{\mu\nu}$ is the Minkowski metric, $H$ is a scalar function encoding the gravitational field strength, and $l_\mu$ is a null vector field. This representation in Cartesian coordinates yields polynomial Christoffel symbols, eliminating the coordinate singularities that plague Boyer-Lindquist representations at the rotation poles.

The Morris-Thorne family covers traversable wormhole geometries through shape function parameterisation. The warp drive family implements Alcubierre-class metrics permitting superluminal effective velocities through spacetime deformation rather than local motion.

## Radiative Transfer

Beyond geodesic tracing, the system implements relativistic radiative transfer for volumetric emission and absorption. The radiative transfer equation along a geodesic is:

$$
\frac{dI_\nu}{d\lambda} = j_\nu - \alpha_\nu I_\nu
$$

where $I_\nu$ is the specific intensity, $j_\nu$ is the emission coefficient, and $\alpha_\nu$ is the absorption coefficient. This enables rendering of accretion disks surrounding black holes with physically-motivated temperature profiles and opacity models.

## Gravitational Lensing

The deflection of light rays near massive objects produces characteristic gravitational lensing phenomena. Strong lensing near black holes creates multiple images and Einstein rings when a background source aligns with the gravitational lens. The photon sphere at $r = 3M$ for Schwarzschild geometry represents the critical surface where light can orbit unstably, producing the characteristic shadow silhouette.

## Architecture

The system is organised into domain-specific components:

| Component | Responsibility |
|-----------|----------------|
| `Sirius.Core` | Domain layer: tensors, metrics, geodesics, physics |
| `Sirius.Render` | GPU acceleration: OptiX, shaders, pipelines |
| `Sirius.Infrastructure` | Application services: UI, configuration |
| `Sirius.Test` | Unit, diagnostic, and benchmark tests |
| `Sirius.Host` | Future development entry points |

## Project Structure

```
Sirius/
├── docs/                    # Comprehensive documentation
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
│   │   ├── OptiX/           # RTX ray tracing
│   │   ├── Shader/          # GLSL shaders
│   │   ├── Pipeline/        # Render pipeline
│   │   ├── Texture/         # Texture assets
│   │   ├── Session/         # Batch rendering
│   │   └── Kernel/          # GPU kernels
│   ├── Sirius.Infrastructure/ # Application services
│   ├── Sirius.Test/         # Test suite
│   └── Sirius.Host/         # Future development
├── bin/                     # Build output
├── ren/                     # Render output (gitignored)
└── CMakeLists.txt           # Build configuration
```

## Documentation

The `docs/` directory contains comprehensive documentation designed to be read in order:

1. **[Philosophy](docs/philosophy.md)**: Design philosophy and worldview
2. **[Foundations](docs/foundations.md)**: Mathematical theory from first principles
3. **[Types](docs/types.md)**: Type system and domain modelling
4. **[Specification](docs/specification.md)**: Formal system requirements
5. **[Guide](docs/guide.md)**: Practical operation
6. **[Standard](docs/standard.md)**: Coding conventions
7. **[Architecture](docs/architecture.md)**: System structure and naming conventions
8. **[Refactor](docs/refactor.md)**: Refactoring discipline
