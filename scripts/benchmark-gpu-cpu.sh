#!/bin/bash
# Benchmark GPU vs CPU Rendering Performance
# Compares rendering times at various resolutions

SIRIUS_BIN="./bin/Sirius.Build/bin/Sirius"
SIRIUS_RENDER="./bin/Sirius.Build/bin/SiriusRender"

if [ ! -f "$SIRIUS_BIN" ]; then
    echo "Error: Sirius binary not found at $SIRIUS_BIN"
    echo "Please build first: cmake --build bin/Sirius.Build"
    exit 1
fi

echo "========================================"
echo "Sirius GPU vs CPU Benchmark"
echo "========================================"
echo ""

# Create output directory
mkdir -p benchmark_output

# Test configurations
declare -a RESOLUTIONS=("256 256" "512 512" "1024 1024")
declare -a SAMPLES=("1" "4" "16")

echo "Testing configurations:"
echo "  Resolutions: 256x256, 512x512, 1024x1024"
echo "  Samples: 1, 4, 16"
echo ""

# Function to run benchmark
run_benchmark() {
    local WIDTH=$1
    local HEIGHT=$2
    local SPP=$3
    local USE_GPU=$4
    local OUTPUT="benchmark_output/bench_${WIDTH}x${HEIGHT}_spp${SPP}"

    if [ "$USE_GPU" = "true" ]; then
        OUTPUT="${OUTPUT}_gpu.ppm"
        GPU_FLAG=""  # GPU is default
    else
        OUTPUT="${OUTPUT}_cpu.ppm"
        GPU_FLAG="--cpu"
    fi

    # Run render and capture time
    local START=$(date +%s.%N)
    $SIRIUS_BIN render --width $WIDTH --height $HEIGHT --spp $SPP --output "$OUTPUT" $GPU_FLAG 2>&1 | grep -E "(Progress|Render|GPU|CPU)" | tail -5
    local END=$(date +%s.%N)

    local ELAPSED=$(echo "$END - $START" | bc)
    echo "  Time: ${ELAPSED}s"
    echo "$ELAPSED"
}

echo "Running benchmarks..."
echo ""

# Results arrays
declare -a GPU_TIMES=()
declare -a CPU_TIMES=()

for RES in "${RESOLUTIONS[@]}"; do
    read -r WIDTH HEIGHT <<< "$RES"

    for SPP in "${SAMPLES[@]}"; do
        echo "----------------------------------------"
        echo "Resolution: ${WIDTH}x${HEIGHT}, Samples: ${SPP}"
        echo ""

        # GPU benchmark (if available)
        echo "GPU Render:"
        GPU_TIME=$(run_benchmark $WIDTH $HEIGHT $SPP "true")
        GPU_TIMES+=("$GPU_TIME")
        echo ""

        # CPU benchmark
        echo "CPU Render:"
        CPU_TIME=$(run_benchmark $WIDTH $HEIGHT $SPP "false")
        CPU_TIMES+=("$CPU_TIME")
        echo ""

        # Calculate speedup
        if [ -n "$GPU_TIME" ] && [ -n "$CPU_TIME" ]; then
            SPEEDUP=$(echo "scale=2; $CPU_TIME / $GPU_TIME" | bc 2>/dev/null || echo "N/A")
            echo "Speedup: ${SPEEDUP}x"
        fi
    done
done

echo ""
echo "========================================"
echo "Benchmark Complete"
echo "========================================"
echo "Output images saved in benchmark_output/"
