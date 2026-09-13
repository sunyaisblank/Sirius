# Recovered retained-precision camera prototypes

These are development sources recovered from recorded tool output from
10 September 2026. They are not called by the production renderer. Keeping
them here preserves the last surviving numerical work in the same Git tree
as Sirius instead of an external test workspace.

`retained_pair.slang` includes the final recorded division-preflight and
square-root range repairs. `recovered-source.json` records the recovered bytes.
The direction module and its two probes precede the complete metric/frame/launch
prototype. The complete launch was still under physical investigation when
the session stopped; its numerical accuracy and performance are unqualified.
Historical reported checks are not current evidence.

The complete launch carries four columns in order: film x/y in pixels and
Cartesian pupil right/up in geometric lengths. Its 104 scientific outputs are
position, tangent, the complete observer frame, four displacements, four
coordinate tangent derivatives, four covariant derivatives, four observer-time
vector derivatives, and a repeated complete observer frame. Every value carries a high part, low part,
absolute error and validity. The complete probe accepts 32 input words and
returns 384 words with metadata, input echoes and 104 retained triples.

Run `python3 tests/support/retained_camera/check_shaders.py` from the checkout
to compile all three probes in both narrow modes. This requires Slang and
SPIR-V Tools. It validates SPIR-V and checks ordered arithmetic, no Float64 or
Int64 arithmetic, denormal preservation and round-to-nearest-even mode. Outputs
stay under `out/retained-camera-build/` and are disposable. This command does
not validate physical execution, error enclosures, camera accuracy, performance,
or integration with transport, events, continuation and detector consumers.

`reference.py` recovers the independent defining metric, general matrix inverse,
observer-frame construction and high-precision differentiation. Its eight
moderate fixtures retain the original calculations; its output paths were moved
under the same ignored directory. With `mpmath==1.3.0` installed in a local
environment, run it to check 100/180-digit stability and the null/frequency
identities and regenerate that reference corpus. This validates the reference
calculation, not shader agreement.

`original_packets.json` preserves the twelve original 68-word scientific input
packets recovered from the September 10 tool output. Run `python3
tests/support/retained_camera/packets.py` to reproduce their binary payload under
`out/retained-camera-build/`. The decoder verifies the original aggregate SHA-256
before writing anything. It preserves float bit patterns; it does not regenerate
the camera samples from decimal approximations. These packets precede the
32-word retained-launch boundary. Pupil coordinates and tangent coefficients
must be observed in the consuming shader context before comparing that launch.

The complete monolithic camera probe currently exceeds 21 GiB of host memory
during pipeline preparation on the pinned WSL2/Dozen Radeon route, including
with Slang `-O0`. Both attempts were stopped before their first dispatch. This
is a preparation failure, with no scientific readback. A bounded staged
implementation is required before further physical qualification; changing
compiler flags alone did not resolve it.
Never discard low parts or error radii to fit the scalar continuation record.
