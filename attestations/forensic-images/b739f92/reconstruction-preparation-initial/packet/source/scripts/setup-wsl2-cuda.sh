#!/bin/bash
# =============================================================================
# WSL2 CUDA Toolkit Setup for Sirius
# =============================================================================
# This script installs CUDA Toolkit in WSL2 Ubuntu.
# Prerequisites: WSL2 with Ubuntu, NVIDIA GPU driver installed on Windows host
# =============================================================================

set -e

echo "=== Sirius WSL2 CUDA Setup ==="
echo ""

# Check if running in WSL2
if [ ! -f /proc/sys/fs/binfmt_misc/WSLInterop ] && [ ! -d /run/WSL ]; then
    echo "ERROR: This script must be run in WSL2"
    exit 1
fi

# Check for NVIDIA driver (exposed through WSL2)
if [ ! -d /usr/lib/wsl/lib ]; then
    echo "ERROR: NVIDIA driver not detected in WSL2"
    echo "Please ensure:"
    echo "  1. NVIDIA GPU driver is installed on Windows host"
    echo "  2. WSL2 GPU support is enabled"
    exit 1
fi

echo "NVIDIA driver detected in WSL2"
nvidia-smi 2>/dev/null || echo "WARNING: nvidia-smi not available yet"

# Add NVIDIA package repository
echo ""
echo "Adding NVIDIA CUDA repository..."
wget -q https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
rm cuda-keyring_1.1-1_all.deb

# Update and install CUDA Toolkit
echo ""
echo "Installing CUDA Toolkit 12.6..."
sudo apt-get update
sudo apt-get install -y cuda-toolkit-12-6

# Set up environment variables
CUDA_PATH="/usr/local/cuda-12.6"
echo ""
echo "Setting up environment variables..."

# Add to .bashrc if not already present
if ! grep -q "CUDA_HOME" ~/.bashrc; then
    cat >> ~/.bashrc << 'EOF'

# CUDA Toolkit
export CUDA_HOME=/usr/local/cuda
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:/usr/lib/wsl/lib:$LD_LIBRARY_PATH
EOF
    echo "Added CUDA environment to ~/.bashrc"
fi

# Source for current session
export CUDA_HOME=/usr/local/cuda
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:/usr/lib/wsl/lib:$LD_LIBRARY_PATH

# Verify installation
echo ""
echo "=== Verification ==="
if command -v nvcc &> /dev/null; then
    echo "CUDA Compiler:"
    nvcc --version
else
    echo "WARNING: nvcc not found. You may need to restart your shell."
fi

echo ""
echo "=== Setup Complete ==="
echo "Please restart your terminal or run: source ~/.bashrc"
echo ""
echo "Next step: Run ./scripts/setup-wsl2-optix.sh to install OptiX SDK"
