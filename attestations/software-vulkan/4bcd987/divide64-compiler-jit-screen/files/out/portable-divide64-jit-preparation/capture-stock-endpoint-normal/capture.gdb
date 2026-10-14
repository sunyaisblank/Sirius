set pagination off
set confirm off
set breakpoint pending on
set debuginfod enabled off
set print thread-events off
set disable-randomization off
set environment LD_PRELOAD /home/astra/.project/Sirius/attestations/software-vulkan/b5d2abf/optimized-ir-stage-code-binding/optimizer_observer.so
set environment SIRIUS_OPT_IR_FOLDER /home/astra/.project/Sirius/out/portable-divide64-jit-preparation/capture-stock-endpoint-normal
source /home/astra/.project/Sirius/out/portable-divide64-jit-preparation/capture_jit.py
source /home/astra/.project/Sirius/out/portable-divide64-jit-preparation/capture_features.py
run
