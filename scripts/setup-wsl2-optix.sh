#!/bin/bash
# =============================================================================
# WSL2 OptiX SDK Setup for Sirius
# =============================================================================
# This script sets up OptiX SDK headers for WSL2 development.
# OptiX is a header-only library; the runtime is provided by the NVIDIA driver.
#
# Prerequisites:
#   - CUDA Toolkit installed (run setup-wsl2-cuda.sh first)
#   - OptiX SDK downloaded from NVIDIA Developer (requires account)
#     https://developer.nvidia.com/designworks/optix/download
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
OPTIX_LOCAL_DIR="$PROJECT_DIR/lib/optix"

echo "=== Sirius WSL2 OptiX Setup ==="
echo ""

# Check if running in WSL2
if [ ! -f /proc/sys/fs/binfmt_misc/WSLInterop ] && [ ! -d /run/WSL ]; then
    echo "ERROR: This script must be run in WSL2"
    exit 1
fi

# Check for CUDA
if ! command -v nvcc &> /dev/null; then
    echo "ERROR: CUDA Toolkit not found"
    echo "Please run ./scripts/setup-wsl2-cuda.sh first"
    exit 1
fi

echo "CUDA found: $(nvcc --version | grep release)"
echo ""

# Check if OptiX is already set up locally
if [ -f "$OPTIX_LOCAL_DIR/include/optix.h" ]; then
    echo "OptiX SDK already installed at: $OPTIX_LOCAL_DIR"
    echo "Version: $(grep '#define OPTIX_VERSION' $OPTIX_LOCAL_DIR/include/optix.h | head -1)"
    exit 0
fi

# Look for OptiX SDK in common locations
OPTIX_SOURCES=(
    "/mnt/c/ProgramData/NVIDIA Corporation/OptiX SDK 9.1.0"
    "/mnt/c/ProgramData/NVIDIA Corporation/OptiX SDK 8.1.0"
    "/mnt/c/ProgramData/NVIDIA Corporation/OptiX SDK 8.0.0"
    "$HOME/NVIDIA-OptiX-SDK-9.1.0-linux64-x86_64"
    "$HOME/NVIDIA-OptiX-SDK-8.1.0-linux64-x86_64"
    "/opt/nvidia/optix"
)

OPTIX_SOURCE=""
for src in "${OPTIX_SOURCES[@]}"; do
    if [ -f "$src/include/optix.h" ]; then
        OPTIX_SOURCE="$src"
        break
    fi
done

if [ -z "$OPTIX_SOURCE" ]; then
    echo "OptiX SDK not found in common locations."
    echo ""
    echo "Please download OptiX SDK from:"
    echo "  https://developer.nvidia.com/designworks/optix/download"
    echo ""
    echo "Options:"
    echo "  1. Install OptiX on Windows host (will be detected at /mnt/c/...)"
    echo "  2. Download Linux version and extract to ~/NVIDIA-OptiX-SDK-*"
    echo "  3. Manually specify path: ./scripts/setup-wsl2-optix.sh /path/to/optix"
    echo ""

    # Check if user provided a path
    if [ -n "$1" ] && [ -f "$1/include/optix.h" ]; then
        OPTIX_SOURCE="$1"
    else
        exit 1
    fi
fi

echo "Found OptiX SDK at: $OPTIX_SOURCE"
echo ""

# Create local optix directory structure
echo "Setting up OptiX headers in project..."
mkdir -p "$OPTIX_LOCAL_DIR"

# Copy only the include directory (headers only, no binaries needed)
cp -r "$OPTIX_SOURCE/include" "$OPTIX_LOCAL_DIR/"

# Verify
if [ -f "$OPTIX_LOCAL_DIR/include/optix.h" ]; then
    VERSION=$(grep '#define OPTIX_VERSION' "$OPTIX_LOCAL_DIR/include/optix.h" | awk '{print $3}')
    MAJOR=$((VERSION / 10000))
    MINOR=$(((VERSION % 10000) / 100))
    PATCH=$((VERSION % 100))

    echo ""
    echo "=== OptiX Setup Complete ==="
    echo "Location: $OPTIX_LOCAL_DIR"
    echo "Version: $MAJOR.$MINOR.$PATCH"
    echo ""
    echo "CMake will automatically detect OptiX at: lib/optix"
else
    echo "ERROR: Failed to set up OptiX headers"
    exit 1
fi
