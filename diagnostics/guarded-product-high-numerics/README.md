# Guarded product-high numerical diagnostics

These are private issue #98 qualification tools, preserved at source revision
4bcd987. They do not install a production arithmetic or shader-policy change.
The candidate header retains host-model globals under `__cplusplus`; it is a
Slang-only input here and must not be used as a production C++ implementation.

Run restored source and execution inputs from
`out/guarded-product-high-numerics`, with staged kernels, recovered executables
and outputs under `bin/linux-gcc/guarded-product-high-numerics`. All original
stock headers other than the private candidate header come from the recorded
Git revision. Exact executed absolute paths, input versions and compiler
commands are retained in the local qualification bundle. Read the input
ledgers before restoring the declared environment; preserve executable modes.

Invoke a declared action with
`python3 out/guarded-product-high-numerics/run_owned.py <action>`.
Every action folder must be new. The runner seals tracked source and declared
inputs, samples owned session memory and host reserve, and closes observed
process births. Do not execute directly from this tracked diagnostics folder,
or reuse an old action folder.

`prepare_residual.py` projects the original ordered residual corpus and four
signed witnesses. The actual shader has no hardware-estimate input.
`compile_numerical.py` keeps the production checker unchanged; its filtered
normal-body text is inspected only, and final unfiltered SPIR-V is assembled,
validated and read again. `verify_products.py` checks typed constants, callee
bodies, fallback, original sum bodies and wire layout.

`verify_residual.py core` and `verify_residual.py rte32` check complete actual
residual words. `verify_numerical.py candidate <mode>` is the original connected
word/interval/refusal/recovery validator. Modes remain `integer`, `normal`,
`normal_transport` and `normal_endpoint`; the private integer sum policy now
contains the guarded FP32 product, with core capabilities retained.

Finite results on the pinned RTE32-capable provider do not test both possible
non-RTE estimates, all platforms, the original eight integrated controls, the
20+2 cost gates, complete science or release qualification. Raw binary evidence
is local and ignored; the evidence tag backs up these sources and recipes.
