#!/bin/bash
# =============================================================================
# WSL2 Development Dependencies for Sirius
# =============================================================================
# Installs build tools and libraries required for Sirius on WSL2 Ubuntu.
# =============================================================================

set -e

echo "=== Sirius WSL2 Dependencies Setup ==="
echo ""

# Check if running in WSL2
if [ ! -f /proc/sys/fs/binfmt_misc/WSLInterop ] && [ ! -d /run/WSL ]; then
    echo "WARNING: Not running in WSL2 (proceeding anyway for native Linux)"
fi

# Update package list
echo "Updating package list..."
sudo apt-get update

# Build essentials
echo ""
echo "Installing build tools..."
sudo apt-get install -y \
    build-essential \
    cmake \
    ninja-build \
    git

# OpenGL and windowing dependencies
echo ""
echo "Installing OpenGL and X11 dependencies..."
sudo apt-get install -y \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libx11-dev \
    libxrandr-dev \
    libxinerama-dev \
    libxcursor-dev \
    libxi-dev \
    libxxf86vm-dev

# Additional utilities
echo ""
echo "Installing additional utilities..."
sudo apt-get install -y \
    pkg-config \
    libwayland-dev \
    libxkbcommon-dev

echo ""
echo "=== Dependencies Installed ==="
echo ""
echo "Next steps:"
echo "  1. ./scripts/setup-wsl2-cuda.sh   - Install CUDA Toolkit"
echo "  2. ./scripts/setup-wsl2-optix.sh  - Set up OptiX SDK"
echo "  3. cmake -B bin/Sirius.Build && cmake --build bin/Sirius.Build"
