#!/usr/bin/env python3
"""Embed the retained SPIR-V admitted for the exact Linux software driver."""

import argparse
import hashlib
from pathlib import Path
import struct


# The ordinary module is part of each selection contract. A renderer with changed
# retained kernels needs a new admission before it can substitute these modules.
STAGES = (
    ("transport", "Transport",
     "2caabcd76c28c5340c21ce5777a58b62b31bb2ed61c1218d15bb10edf0304752",
     "1e48f26028368c7b6b20ec5e7fe0b7283ee084327e7a9d6f22b587ec6e3dda4d"),
    ("endpoint", "Endpoint",
     "d7c01a3d67f7c4617643aa36c2369ffc190f2e36e2ec9fd1d1aeb0a346a45e7c",
     "e505643b3ec8865702d035f062139a07da73a8e485210fd7ca094e20d2a621c5"),
    ("dopri_phase", "DopriPhase",
     "10d94ce407e6ad156920e9cde616f74c84ab88606377e5151e4115c8ce08e4a3",
     "f26d03f631a22e025b9a20ba700e86ef3093630935817e4611c3552094c1a336"),
)


def checked_spirv(path, expected):
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != expected:
        raise ValueError(f"{path} does not match the qualified retained module SHA-256")
    if len(payload) % 4 or payload[:4] != b"\x03\x02\x23\x07":
        raise ValueError(f"{path} is not little-endian SPIR-V")
    return payload


def emit(args):
    lines = [
        '#pragma once',
        '#include "retained_kernels.h"',
        'namespace sirius::backend::retained_program {',
        f'inline constexpr char kPackagedDriverSha256[] = "{args.driver_sha256}";',
        f'inline constexpr char kPackagedLlvmSha256[] = "{args.llvm_sha256}";',
    ]
    for stage, name, ordinary, sealed in STAGES:
        checked_spirv(args.ordinary / f"retained_{stage}.spv", ordinary)
        payload = checked_spirv(args.bundle / f"{stage}.spv", sealed)
        words = [word for word, in struct.iter_unpack("<I", payload)]
        lines.append(f'inline constexpr char k{name}SealedSha256[] = "{sealed}";')
        lines.append(f'inline constexpr std::array<std::uint32_t,{len(words)}> k{name}SealedShader{{{{')
        for index in range(0, len(words), 12):
            lines.append(','.join(f'{word}u' for word in words[index:index + 12]) + ',')
        lines.append('}};')
    lines.append('}')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pending = args.output.with_suffix('.pending')
    pending.write_text('\n'.join(lines) + '\n', encoding='ascii')
    pending.replace(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--ordinary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--driver-sha256', required=True)
    parser.add_argument('--llvm-sha256', required=True)
    try:
        emit(parser.parse_args())
    except (OSError, ValueError) as error:
        parser.exit(1, f"{error}\n")
