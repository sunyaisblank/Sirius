Completed at clean pushed source [6d2af336](https://github.com/sunyaisblank/Sirius/commit/6d2af336b12dc430792bd5a83c04fb5632868c04). Assigning the complete validated CameraState clears all six held movement flags on successful reinitialisation. The original input control now checks each held key across Stop/Initialise, unchanged idle coordinates/no restart, and fresh input recovery; both strict GCC14 and Clang21 runs pass one case with zero failures/errors/skips/disabled.

Both existing strict C++26 presets build all seven roles and source/alignment gates. All 35 retained SPVs and three generated math headers remain byte-exact. CI39487479545 passes all three platform integration jobs/30 authority controls; the draft PR's CI39487799597 also passes, with actual merge tree equal to source6d2af336. Five full-qualification/diagnostic jobs per run were unrequested.

[Complete preserved source/counterexample/products/owner/CI records](https://github.com/sunyaisblank/Sirius/tree/ad36d3bd6d8ffabcef426c5d8429e3b617f9b45d/attestations/qualification/6d2af336/viewer-input-reinitialization). This closes the input-reset defect only. No native input-delivery, frame-speed, full Mandatory/scientific/provider/release qualification follows; current external-domain admission remains0/8. PR62 stays draft, main stays5e706073 and #44/#48/#53 remain open.

Original finding and acceptance retained:

`InteractiveViewer::Initialise` overwrites only the four camera coordinates. Six held movement flags survive a stopped viewer, so a missing key-release event in the previous window can move a freshly initialised camera before any new input. This is an established source/public-API counterexample, not a native-event or frame-speed claim.

The exact existing GCC strict CMake app library at 9b770714 reproduces `held_motion_after_reinitialization=1` after pressing all six movement keys, Stop and successful Initialise. The bounded source-only probe completed in 5.401 seconds with no provider/render launch.

Repair: assign the already validated complete initial CameraState, which includes false held-key flags, under the existing state mutex. Extend the existing press/repeat/release/mouse/scroll control with stopped reinitialisation for each movement key, stable idle coordinates/no restart, and fresh-input recovery. Retain all original assertions, source guards, physics, dimensions, budgets and registration. Invalid configuration continues to return before state changes.

Validate through existing strict GCC/Clang presets and the existing non-render app control. Preserve the source/result, identify the development branch, and keep full native window/input, Mandatory, scientific/provider/frame/main/install/release acceptance separate. #44/#48/#53 and draft PR62 remain open.
