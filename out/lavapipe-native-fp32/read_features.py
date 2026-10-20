"""Read original raw argument captures; never creates a device or alters old evidence."""
from pathlib import Path
import hashlib
import json
import struct


def inspect(folder):
    rows = json.loads((folder / "device-features.json").read_text())
    assert len(rows["device_creations"]) == len(rows["shader_modules"]) == 1
    device, shader = rows["device_creations"][0], rows["shader_modules"][0]
    assert device["completed"] and device["result"] == 0 and device["handle"]
    assert shader["completed"] and shader["result"] == 0 and shader["handle"]
    assert shader["device"] == device["handle"]
    actual = device["features"]
    raw = bytes.fromhex(actual["raw_hex"])
    assert len(raw) == 220 and hashlib.sha256(raw).hexdigest() == actual["sha256"]
    assert list(struct.unpack("<55I", raw)) == list(actual["values"].values())
    values = actual["values"]
    assert values["shaderInt64"] == 1
    assert all(value == 0 for name, value in values.items()
               if name not in ("shaderFloat64", "shaderInt64"))
    queries = [q for q in rows["queries"] if q["physical_device"] == device["physical_device"]]
    assert queries and all(q["completed"] for q in queries)
    for query in queries:
        qraw = bytes.fromhex(query["features"]["raw_hex"])
        assert len(qraw) == 220 and hashlib.sha256(qraw).hexdigest() == query["features"]["sha256"]
        assert list(struct.unpack("<55I", qraw)) == list(query["features"]["values"].values())
        assert query["features"]["values"]["shaderInt64"] == 1
        assert query["features"]["values"]["shaderFloat64"] == values["shaderFloat64"]
    loaded = json.loads((folder / "actual-loaded-module.json").read_text())
    assert shader["code_sha256"] == loaded["sha256"] and shader["code_bytes"] == loaded["bytes"]
    assert device["flags"] == 0 and device["queue_count"] == 1 and device["enabled_layer_count"] == 0
    assert device["queue_priorities_count"] == 1 and device["queue_priority"] == 1.0
    for item in device["feature_chain"]:
        assert item["sType"] == 1000579000 and item["fma16_32_64"] == [0, 1, 0]
    return {"actual_public_create_arguments_enable_shaderInt64": True,
            "successful_logical_device_used_for_actual_complete_shader": True,
            "features": values, "queue_family": device["queue_family"],
            "enabled_extension_count": device["enabled_extension_count"],
            "enabled_extensions": device["enabled_extensions"],
            "feature_chain_values": [item["fma16_32_64"] for item in device["feature_chain"]],
            "scope": "Observed public loader arguments and successful device/shader creation only; no active retained arithmetic or benefit verdict."}
