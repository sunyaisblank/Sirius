set pagination off
set confirm off
set breakpoint pending on
set debuginfod enabled off
set print thread-events off
set disable-randomization off
set environment LD_PRELOAD /home/astra/.project/Sirius/attestations/software-vulkan/b5d2abf/optimized-ir-stage-code-binding/optimizer_observer.so
set environment SIRIUS_OPT_IR_FOLDER /home/astra/.project/Sirius/out/portable-divide-fp32-estimate-restoring-jit/capture-candidate-endpoint-normal
source /home/astra/.project/Sirius/out/portable-divide-fp32-estimate-restoring-jit/capture_jit.py
source /home/astra/.project/Sirius/out/portable-divide-fp32-estimate-restoring-jit/capture_features.py
run
