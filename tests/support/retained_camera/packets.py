"""Decode the recovered scientific packets without decimal float conversion."""

import hashlib
import json
from pathlib import Path
import struct


SOURCE = Path(__file__).resolve().parent
OUTPUT = SOURCE.parents[2] / "out" / "retained-camera-build"
EXPECTED_SHA256 = "f4c7735c37a0030adedc1aeeae862a49787fcc5bfbd9e9c1322724f482231c4f"


def original_packets():
    document = json.loads((SOURCE / "original_packets.json").read_text())
    if document["schema"] != "sirius-original-camera-packets-v1":
        raise ValueError("unknown packet schema")
    base = [int(word, 16) for word in document["base_float32_bits"]]
    indices = document["variable_indices"]
    if len(base) != 68 or indices != [31, 32, 44, 45, 46, 66, 67]:
        raise ValueError("scientific packet layout changed")
    names, packets = [], []
    for case in document["cases"]:
        words = base.copy()
        values = case["variable_float32_bits"]
        if len(values) != len(indices):
            raise ValueError("incomplete scientific case")
        for index, word in zip(indices, values):
            words[index] = int(word, 16)
        names.append(case["name"])
        packets.append(struct.pack("<68I", *words))
    if len(packets) != 12 or len(set(names)) != 12:
        raise ValueError("expected twelve distinct original cases")
    digest = hashlib.sha256(b"".join(packets)).hexdigest()
    if digest != EXPECTED_SHA256 or digest != document["aggregate_sha256"]:
        raise ValueError("scientific bytes differ from the historical manifest")
    return names, packets


def main():
    names, packets = original_packets()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "original-scientific-packets.bin").write_bytes(b"".join(packets))
    print(f"Recovered {len(names)} exact 68-word packets; SHA-256 {EXPECTED_SHA256}")


if __name__ == "__main__":
    main()
