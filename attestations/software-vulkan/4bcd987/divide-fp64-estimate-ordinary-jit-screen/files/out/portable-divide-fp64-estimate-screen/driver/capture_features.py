"""Passive addon to the unchanged x86_64 GDB capture; record public Vulkan arguments."""
import gdb
import hashlib
import json
import os
import pathlib
import struct

feature_folder = pathlib.Path(os.environ["SIRIUS_GEOMETRY_CAPTURE"])
feature_layout = json.loads(pathlib.Path(os.environ["SIRIUS_FEATURE_LAYOUT"]).read_text())
feature_loader = pathlib.Path(os.environ["SIRIUS_FEATURE_LOADER"]).resolve(strict=True)
feature_rows = {"queries": [], "device_creations": [], "shader_modules": []}


def feature_write():
    (feature_folder / "device-features.json").write_text(json.dumps(feature_rows, indent=2) + "\n")


def feature_boundary():
    pc = reg("rip")
    library = gdb.solib_name(pc)
    assert library and pathlib.Path(library).resolve(strict=True) == feature_loader, (pc, library)
    return {"pc": pc, "symbol": gdb.newest_frame().name(), "library": library,
            "maps": pathlib.Path("/proc/" + str(gdb.selected_inferior().pid) + "/maps").read_text()}


def feature_bytes(pointer):
    raw = memory(pointer, feature_layout["features_bytes"])
    values = struct.unpack("<55I", raw)
    return {"pointer": pointer, "raw_hex": raw.hex(), "sha256": hashlib.sha256(raw).hexdigest(),
            "values": dict(zip(feature_layout["features_fields"], values))}


def feature_string(pointer):
    value = bytearray()
    for offset in range(256):
        byte = memory(pointer + offset, 1)
        if byte == b"\0":
            return value.decode("ascii")
        value.extend(byte)
    raise ValueError("unterminated Vulkan extension name")


class FeaturesReturned(gdb.FinishBreakpoint):
    def __init__(self, row, pointer):
        super().__init__(gdb.newest_frame(), internal=True)
        self.row, self.pointer = row, pointer

    def stop(self):
        try:
            self.row.update(completed=True, features=feature_bytes(self.pointer))
            feature_write()
            return False
        except Exception as error:
            return failure(error)


class FeaturesQuery(gdb.Breakpoint):
    def stop(self):
        try:
            row = {"physical_device": reg("rdi"), "boundary": feature_boundary(), "completed": False}
            feature_rows["queries"].append(row)
            FeaturesReturned(row, reg("rsi"))
            feature_write()
            return False
        except Exception as error:
            return failure(error)


class HandleReturned(gdb.FinishBreakpoint):
    def __init__(self, row, pointer):
        super().__init__(gdb.newest_frame(), internal=True)
        self.row, self.pointer = row, pointer

    def stop(self):
        try:
            result = struct.unpack("<i", struct.pack("<I", reg("rax") & 0xffffffff))[0]
            self.row.update(completed=True, result=result,
                            handle=struct.unpack("<Q", memory(self.pointer, 8))[0] if result == 0 else None)
            feature_write()
            return False
        except Exception as error:
            return failure(error)


class DeviceCreation(gdb.Breakpoint):
    def stop(self):
        try:
            pointer = reg("rsi")
            raw = memory(pointer, 72)
            u32 = lambda at: struct.unpack_from("<I", raw, at)[0]
            u64 = lambda at: struct.unpack_from("<Q", raw, at)[0]
            assert u32(0) == 3 and u32(20) == 1 and u32(32) == 0
            chain, next_pointer = [], u64(8)
            while next_pointer:
                assert len(chain) < 4 and next_pointer not in [item["pointer"] for item in chain]
                header = memory(next_pointer, 16)
                kind, following = struct.unpack_from("<I", header)[0], struct.unpack_from("<Q", header, 8)[0]
                assert kind == 1000579000, "unexpected feature chain in this exact source route"
                value = memory(next_pointer, 32)
                chain.append({"pointer": next_pointer, "sType": kind, "raw_hex": value.hex(),
                              "fma16_32_64": list(struct.unpack_from("<3I", value, 16))})
                next_pointer = following
            queue = memory(u64(24), 40)
            queue_priority = struct.unpack("<f", memory(struct.unpack_from("<Q", queue, 32)[0], 4))[0]
            assert u32(48) < 32
            extensions = [feature_string(struct.unpack("<Q", memory(u64(56) + 8 * i, 8))[0])
                          for i in range(u32(48))]
            row = {"physical_device": reg("rdi"), "boundary": feature_boundary(),
                   "create_info_pointer": pointer, "raw_create_info_hex": raw.hex(),
                   "sType": u32(0), "flags": u32(16), "queue_count": u32(20),
                   "enabled_layer_count": u32(32), "enabled_extension_count": u32(48),
                   "enabled_extensions": extensions,
                   "raw_queue_info_hex": queue.hex(), "queue_family": struct.unpack_from("<I", queue, 20)[0],
                   "queue_priorities_count": struct.unpack_from("<I", queue, 24)[0],
                   "queue_priority": queue_priority, "feature_chain": chain,
                   "features": feature_bytes(u64(64)), "completed": False}
            feature_rows["device_creations"].append(row)
            HandleReturned(row, reg("rcx"))
            feature_write()
            return False
        except Exception as error:
            return failure(error)


class ShaderCreation(gdb.Breakpoint):
    def stop(self):
        try:
            raw = memory(reg("rsi"), 40)
            size, pointer = struct.unpack_from("<QQ", raw, 24)
            assert 20 < size < 8 * 1024**2 and size % 4 == 0
            code = memory(pointer, size)
            assert code == pathlib.Path(os.environ["SIRIUS_MODULE_INPUT"]).read_bytes()
            row = {"device": reg("rdi"), "boundary": feature_boundary(),
                   "raw_create_info_hex": raw.hex(), "code_bytes": size,
                   "code_sha256": hashlib.sha256(code).hexdigest(), "completed": False}
            feature_rows["shader_modules"].append(row)
            HandleReturned(row, reg("rcx"))
            feature_write()
            return False
        except Exception as error:
            return failure(error)


FeaturesQuery("vkGetPhysicalDeviceFeatures", internal=True)
DeviceCreation("vkCreateDevice", internal=True)
ShaderCreation("vkCreateShaderModule", internal=True)
feature_write()
