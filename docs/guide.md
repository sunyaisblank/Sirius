# Practical Guide

*From Configuration to Rendering*

> "In theory, there is no difference between theory and practice. In practice, there is."
> — Attributed to Yogi Berra

## Preface

This document bridges the gap between theory and operation. It covers practical patterns for building, configuring, and running the Sirius renderer, troubleshooting common issues, and optimising performance.

The reader is assumed to have read the [Philosophy](philosophy.md) and [Foundations](foundations.md) documents, to understand basic general relativity concepts, and to be comfortable with command-line interfaces.

---

## Part I: Building Sirius

### 1.1 Prerequisites

| Requirement | Version | Notes |
|-------------|---------|-------|
| CMake | ≥ 3.20 | Build system |
| C++ Compiler | C++17 support | GCC 9+, Clang 10+, MSVC 2019+ |
| CUDA Toolkit | ≥ 11.0 | For GPU acceleration |
| OptiX | 7.x | NVIDIA ray tracing |
| OpenGL | 4.5+ | Display |
| GLFW | 3.3+ | Window management |
| GLM | 0.9.9+ | Maths library |
| Dear ImGui | 1.89+ | UI |

### 1.2 Build Configuration

```bash
# Clone repository
git clone https://github.com/username/sirius.git
cd sirius

# Create build directory
mkdir Sirius.Build
cd Sirius.Build

# Configure
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build . --config Release
```

### 1.3 Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `SIRIUS_ENABLE_TESTS` | ON | Build test suite |
| `SIRIUS_ENABLE_BENCHMARKS` | ON | Build benchmarks |
| `SIRIUS_CUDA_ARCH` | Auto | Target GPU architecture |
| `SIRIUS_DOUBLE_PRECISION` | OFF | Use double precision on GPU |

---

## Part II: Running the Renderer

### 2.1 Basic Usage

```bash
# Interactive mode (default)
./SiriusRender

# With configuration file
./SiriusRender --config config.json

# Offline render
./SiriusOffline --config render_config.json --output image.exr
```

### 2.2 Command-Line Options

| Option | Description |
|--------|-------------|
| `--config <file>` | Load configuration from JSON file |
| `--metric <name>` | Select metric (schwarzschild, kerr, etc.) |
| `--resolution <WxH>` | Set resolution (e.g., 1920x1080) |
| `--output <file>` | Output file path (offline mode) |
| `--help` | Show help message |

### 2.3 Interactive Controls

| Key | Action |
|-----|--------|
| WASD | Move camera |
| Mouse | Rotate view |
| Scroll | Adjust distance |
| Space | Pause/resume |
| F1 | Toggle UI |
| F5 | Reload shaders |
| F12 | Screenshot |

---

## Part III: Metric Configuration

### 3.1 Available Metrics

| Metric | Code | Parameters |
|--------|------|------------|
| Minkowski | `minkowski` | None |
| Schwarzschild | `schwarzschild` | Mass |
| Kerr | `kerr` | Mass, Spin |
| Reissner-Nordström | `reissner_nordstrom` | Mass, Charge |
| Kerr-Newman | `kerr_newman` | Mass, Spin, Charge |
| Ellis Wormhole | `ellis` | Throat radius |
| Alcubierre | `alcubierre` | Velocity, Radius, Thickness |

### 3.2 Schwarzschild Configuration

```json
{
  "metric": {
    "type": "schwarzschild",
    "mass": 1.0
  }
}
```

Key properties:

- Event horizon at $r = 2M$
- Photon sphere at $r = 3M$
- ISCO at $r = 6M$

### 3.3 Kerr Configuration

```json
{
  "metric": {
    "type": "kerr",
    "mass": 1.0,
    "spin": 0.9
  }
}
```

The spin parameter $a$ is specified as a fraction of maximum ($a_{max} = M$). Valid range: $[0, 0.998]$.

Effects of increasing spin:

- Horizons move closer together
- Ergosphere expands
- Frame dragging intensifies
- Shadow becomes asymmetric

---

## Part IV: Camera Configuration

### 4.1 Observer Setup

```json
{
  "camera": {
    "distance": 20.0,
    "inclination": 85.0,
    "azimuth": 0.0,
    "fov": 60.0
  }
}
```

| Parameter | Description | Range |
|-----------|-------------|-------|
| `distance` | Distance from origin in units of M | [5, 1000] |
| `inclination` | Angle from polar axis in degrees | [0, 180] |
| `azimuth` | Angle around axis in degrees | [0, 360] |
| `fov` | Field of view in degrees | [10, 120] |

### 4.2 Views for Common Effects

| Effect | Inclination | Distance |
|--------|-------------|----------|
| Face-on disk | 0° | 20M |
| Edge-on disk | 90° | 30M |
| Shadow emphasis | 85° | 15M |
| Einstein ring | Variable | Aligned with source |

---

## Part V: Accretion Disk

### 5.1 Disk Configuration

```json
{
  "disk": {
    "enabled": true,
    "innerRadius": 6.0,
    "outerRadius": 20.0,
    "temperatureScale": 1.0,
    "emissionModel": "novikov_thorne"
  }
}
```

### 5.2 Emission Models

| Model | Description |
|-------|-------------|
| `uniform` | Constant emission |
| `power_law` | $T \propto r^{-3/4}$ |
| `novikov_thorne` | Relativistic thin disk |

### 5.3 Physical Parameters

| Parameter | Description |
|-----------|-------------|
| `innerRadius` | Inner edge (usually ISCO) |
| `outerRadius` | Outer extent of disk |
| `temperatureScale` | Multiplier for temperature profile |
| `Mdot` | Accretion rate (for absolute temperatures) |

---

## Part VI: Performance Tuning

### 6.1 Quality vs Speed

| Setting | Performance Impact | Quality Impact |
|---------|-------------------|----------------|
| Resolution | Linear | Direct |
| Max steps | Linear | Accuracy near horizon |
| Step size | Inverse | Integration accuracy |
| Adaptive stepping | Slight overhead | Significant accuracy gain |

### 6.2 Recommended Presets

**Fast Preview:**

```json
{
  "render": {
    "resolution": [960, 540],
    "maxSteps": 500,
    "adaptiveStepping": false
  }
}
```

**Quality Render:**

```json
{
  "render": {
    "resolution": [1920, 1080],
    "maxSteps": 2000,
    "adaptiveStepping": true,
    "tolerance": 1e-8
  }
}
```

**Publication Quality:**

```json
{
  "render": {
    "resolution": [3840, 2160],
    "maxSteps": 5000,
    "adaptiveStepping": true,
    "tolerance": 1e-10,
    "samples": 16
  }
}
```

### 6.3 GPU Memory Management

| Resolution | Approximate VRAM |
|------------|------------------|
| 1080p | ~500 MB |
| 1440p | ~750 MB |
| 4K | ~1.5 GB |

---

## Part VII: Troubleshooting

### 7.1 Common Issues

**Black screen output**

Possible causes:

- Observer inside event horizon
- Rays terminated before escaping
- Background texture not loaded

Resolution:

- Increase observer distance
- Increase max integration steps
- Check texture path in configuration

**Visible artefacts**

Possible causes:

- Insufficient integration steps
- Coordinate singularity intersection
- Floating-point precision limits

Resolution:

- Increase `maxSteps`
- Enable adaptive stepping
- Check ray termination conditions

**Low frame rate**

Possible causes:

- High resolution
- Complex metric (numerical)
- Excessive integration steps

Resolution:

- Reduce resolution or render scale
- Use analytic metric if applicable
- Tune step size and max steps

### 7.2 Diagnostic Mode

```bash
./SiriusRender --diagnostic
```

Outputs:

- Conservation law violations per frame
- NaN/Inf ray count
- Mean steps per ray
- GPU memory usage

### 7.3 Log Analysis

Logs are written to `sirius.log`:

```
[INFO] Frame 100: 2073600 rays, 847 steps/ray avg, 45.2 FPS
[WARN] Ray (512, 384): Conservation error 0.0012 exceeds tolerance
[ERROR] Ray (100, 200): NaN detected in position
```

---

## Part VIII: Offline Rendering

### 8.1 Batch Rendering

```bash
./SiriusOffline --config render.json --output frame.exr
```

### 8.2 Animation Sequences

```bash
# Render sequence with varying spin
for spin in $(seq 0 0.1 0.99); do
  ./SiriusOffline --metric kerr --spin $spin --output "frame_${spin}.exr"
done
```

### 8.3 Output Formats

| Format | Extension | Use Case |
|--------|-----------|----------|
| OpenEXR | `.exr` | HDR, post-processing |
| PNG | `.png` | Web, presentation |
| PPM | `.ppm` | Simple, portable |

---

## Appendix A: Configuration Reference

### Complete Configuration Example

```json
{
  "metric": {
    "type": "kerr",
    "mass": 1.0,
    "spin": 0.9
  },
  "camera": {
    "distance": 20.0,
    "inclination": 85.0,
    "azimuth": 0.0,
    "fov": 60.0
  },
  "disk": {
    "enabled": true,
    "innerRadius": 6.0,
    "outerRadius": 20.0,
    "emissionModel": "novikov_thorne"
  },
  "render": {
    "resolution": [1920, 1080],
    "maxSteps": 1500,
    "stepSize": 0.1,
    "adaptiveStepping": true,
    "tolerance": 1e-6
  },
  "background": {
    "texture": "starfield.png"
  }
}
```

---

*End of Guide*
