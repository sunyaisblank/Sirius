# Saved SQTT detailed-target boundary at 6be3132 (#44)

Transport's 436-byte allocated tail is exactly 109 GFX11 `s_code_end` words. Binding those proven boundaries resolves the requested entry-end address, but the 0.1.3 wave arrays remain byte-identical and auxiliary address `0x10001b800` remains unbound.

The saved detailed target is WGP 0, SIMD 0. Only the first Transport and an auxiliary wave start on that target; no Endpoint wave does. Expensive later Transport waves occupy other lanes. The pinned RDP API/CLI/serializer and PAL GpaSession path expose no supported lane selector. Lower PAL options do not establish a supported per-owned-RDP control. No guessed keys or driver settings are accepted.

An actual bounded CPU-only official decoder 0.1.6 call recovers 2,663 anonymous memory-operation records on WGP 0: 1,077 VMEM and 1,586 LDS. They have no PC, wave or SIMD identity, and execution cycles exclude stall. They cannot establish exclusive source/event cost. Both returned wave timeline/end/occupancy intervals now agree, while the auxiliary wave still has 16 zero PCs. No INFO record is not complete source coverage.

The same 11,912 ordered reference-clock values survive, but every shader-clock value changes under the newer decoder. Shader-array encoding also changes. Keep both versions distinct; the earlier SPM clock comparison is historical and cannot be reused for this result. Neither version yields a reference-clock frequency or exact SPM/event/display origin.

No optimization is adopted. Original source/products, complete numerical and ordered-stream contracts, benefit gates and rendering requirements remain authoritative. Useful source and receipts are preserved in pushed Git; detailed raw evidence remains local. No new GPU dispatch, rendered frame, scientific or release qualification was run in this batch. #44, #48 and #53 remain open.
