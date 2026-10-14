#!/bin/bash
# Render test script for Sirius black holes
# Run this from the Sirius project directory

SIRIUS="./bin/Sirius.Build/bin/Sirius"

if [ ! -f "$SIRIUS" ]; then
    echo "Error: Sirius binary not found at $SIRIUS"
    echo "Please build first: cmake --build bin/Sirius.Build"
    exit 1
fi

echo "=== Sirius Black Hole Render Test ==="
echo ""

# Create output directory
mkdir -p renders

# Render Schwarzschild (non-spinning black hole)
echo "Rendering Schwarzschild black hole..."
$SIRIUS render \
    -m Schwarzschild \
    -w 512 -h 512 \
    -s 32 \
    -d 30 \
    -i 80 \
    --fov 60 \
    -o renders/schwarzschild.png

if [ $? -eq 0 ]; then
    echo "  -> renders/schwarzschild.png created"
else
    echo "  -> Failed to render Schwarzschild"
fi

echo ""

# Render Kerr (spinning black hole, a=0.9)
echo "Rendering Kerr black hole (spin=0.9)..."
$SIRIUS render \
    -m Kerr \
    -a 0.9 \
    -w 512 -h 512 \
    -s 32 \
    -d 30 \
    -i 80 \
    --fov 60 \
    -o renders/kerr.png

if [ $? -eq 0 ]; then
    echo "  -> renders/kerr.png created"
else
    echo "  -> Failed to render Kerr"
fi

echo ""
echo "=== Done ==="
echo "Check the 'renders' directory for output images."
