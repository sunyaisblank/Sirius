# Rebuild baseline tapes

Reference renders for the byte-identity gate. Untracked; regenerate with the
commands below and compare with `cmp` after any change that claims to leave
the render path untouched.

Provenance: linux-gcc preset (GCC 14, Release), regenerated 2026-07-18 after
the tracer step-rejection fix (a rejected RK45 step now retries at the reduced
step instead of terminating the ray as MaxSteps). That fix moved 16.1% of
pixels in the Kerr scene, all inside the strong-field bounding box around the
shadow and photon ring (rows 141-367, cols 124-376 of 512^2; the far field was
byte-identical); the Alcubierre scene was byte-identical pre/post. The prior
tapes (640x360 kerr, 320x180 warp, taped 2026-07-17 at commit 58c97da) had no
reproducing script and are superseded by these documented geometries.

These are CPU-reference tapes; since backend auto go-live (2026-07-18)
resolves to Vulkan when a device is visible, every command pins `--cpu`.

    SIRIUS=./bin/linux-gcc/src/sirius/app/sirius
    $SIRIUS render --cpu -m Kerr -a 0.9 -w 512 -h 512 -s 32 -d 30 -i 80 --fov 60 \
        -o renders/rebuild-baseline/kerr.png     # and kerr.exr
    $SIRIUS render --cpu -m Alcubierre --warp-velocity 0.5 --bubble-radius 1.0 \
        -w 512 -h 512 -s 16 -d 20 -i 80 --fov 90 \
        -o renders/rebuild-baseline/warp.png
    $SIRIUS render --cpu -m Morris-Thorne --throat-radius 1.0 \
        -w 512 -h 512 -s 16 -d 15 -i 80 --fov 60 \
        -o renders/rebuild-baseline/wormhole.png
