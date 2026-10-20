# Viewer input reinitialisation

Successful `Initialise` now assigns the complete validated camera state, clearing all six held movement keys. A release event lost when a prior window stops cannot steer a freshly initialised viewer. Invalid configuration still returns before state changes. Existing running-state admission, locking, physics and rendering settings remain unchanged.

At canonical baseline `9b770714`, an owned public-API counterexample linked to the existing strict app library prints `held_motion_after_reinitialization=1`. It completed in 5.401 seconds; this confirms stale input state, not native event delivery.

At clean pushed source `6d2af336`, GCC14 and Clang21 strict C++26 qualification presets build all seven executable roles and source/alignment gates. Incremental build owners pass in 38.207 and 30.912 seconds, with 1,038,116 and 637,648 KiB sampled aggregate peaks, original 1800-second/4-GiB guards, 2-GiB reserve and two build jobs. The existing press/repeat/release/mouse/scroll control additionally checks every movement key across stopped reinitialisation, stable idle coordinates/no restart, and fresh input recovery. Both executions pass one case with zero failures/errors/skips/disabled cases. The short controls ended before a nonzero RSS sample; no zero-memory claim is made.

Both presets retain 1,152 registered cases and 0/8 external-domain admission. Registration is not execution. All 35 retained SPVs and three generated math headers equal the preserved canonical 9b baseline byte-for-byte on both compilers. This change justifies no repeated physics/benefit/failed-frame study.

The deduplicated archive preserves all 234 logical payload paths, with exact byte/mode backings for 138 unique members: source/patch, original baseline probe and selectors, bounded owner/input/log/XML records, actual local products, registration, copied resource volumes and exact math comparison. Builds are incremental; unchanged dependencies retain earlier producer bytes. Platform CI terminal results are separate. CI artifact API/container identities alone do not back up their native executable contents.

The full retained frame, native window/input, scientific, Mandatory/sanitizer, external-domain and exact-main installed-release requirements remain open. PR62 remains draft and main stays 5e706073.

CI 39487479545 and PR CI 39487799597 both finish SUCCESS: three strict platform integration jobs and 30 original authority controls per run. Five full-qualification/diagnostic jobs per run were unrequested. The actual PR merge fc3d5ec8 has the complete 6d2af336 source tree and retains its own product revision. Complete job/API/log identities are in `ci-review.json`; native artifact container/API metadata is recorded without claiming a local backup or admission of its contents.
