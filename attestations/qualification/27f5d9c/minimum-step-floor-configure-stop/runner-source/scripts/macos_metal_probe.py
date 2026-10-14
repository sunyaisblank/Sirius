#!/usr/bin/env python3
"""One bounded Initialize MSL library/function/PSO diagnostic; no dispatch.

Run only under the existing Darwin owner (300 s / 4 GiB RSS / 2 GiB reserve).
Uses synchronous error-returning Metal methods instead of MoltenVK's asynchronous
callbacks, with source-derived compile options and compute descriptor. Optional
Size changes only the library optimization setting as a diagnostic experiment.
An explicit source control removes only the captured forced-inline annotations
while retaining the selected options and every function body; both files stay sealed.
No product qualification or equivalence of the two compilation routes is claimed.
"""

import argparse
import ctypes as C
import hashlib
import json
from pathlib import Path
import platform
import sys
import time


MSL_SHA = "b2d7daec9928356de4ca24a724d135adc59685a330b043e4fa55f2fe2f9b5550"
MSL_BYTES = 85091
SOURCE = "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc"
FORCED_INLINE = b"static inline __attribute__((always_inline))"
FORCED_INLINE_COUNT = 84


def require(ok, message):
    if not ok:
        raise ValueError(message)


def unforced_msl(original):
    """Derive only this frozen source control, never accept arbitrary MSL."""
    require(len(original) == MSL_BYTES and hashlib.sha256(original).hexdigest() == MSL_SHA,
            "MSL differs from exact captured Initialize source")
    require(original.count(FORCED_INLINE) == FORCED_INLINE_COUNT
            and original.count(b"__attribute__((always_inline))") == FORCED_INLINE_COUNT,
            "captured forced-inline declaration inventory differs")
    return original.replace(FORCED_INLINE, b"static inline")


def emit(event, **values):
    print(json.dumps({"event": event, "monotonic_seconds": time.monotonic(),
                      "utc_ns": time.time_ns(), **values}), flush=True)


class ObjC:
    """Typed objc_msgSend calls; no structs, floating-point ABI or callbacks."""

    def __init__(self):
        self.runtime = C.CDLL("/usr/lib/libobjc.A.dylib")
        self.foundation = C.CDLL("/System/Library/Frameworks/Foundation.framework/Foundation")
        self.metal = C.CDLL("/System/Library/Frameworks/Metal.framework/Metal")
        for name in ["objc_getClass", "sel_registerName"]:
            fn = getattr(self.runtime, name)
            fn.argtypes, fn.restype = [C.c_char_p], C.c_void_p
        self.runtime.object_getClassName.argtypes = [C.c_void_p]
        self.runtime.object_getClassName.restype = C.c_char_p
        self.address = C.cast(self.runtime.objc_msgSend, C.c_void_p).value
        self.calls = {}

    def cls(self, name):
        result = self.runtime.objc_getClass(name.encode())
        require(result, "Objective-C class missing: " + name)
        return result

    def send(self, obj, selector, result=C.c_void_p, types=(), args=()):
        require(obj, "nil receiver: " + selector)
        key = (result, tuple(types))
        if key not in self.calls:
            self.calls[key] = C.CFUNCTYPE(result, C.c_void_p, C.c_void_p, *types)(self.address)
        return self.calls[key](obj, self.runtime.sel_registerName(selector.encode()), *args)

    def integer(self, obj, selector):
        return self.send(obj, selector, C.c_size_t)

    def boolean(self, obj, selector):
        return bool(self.send(obj, selector, C.c_bool))

    def supports(self, obj, selector):
        return bool(self.send(obj, "respondsToSelector:", C.c_bool, (C.c_void_p,),
                              (self.runtime.sel_registerName(selector.encode()),)))

    def empty_collection(self, obj, selector):
        value = self.send(obj, selector)
        return not value or self.integer(value, "count") == 0

    def string(self, text):
        raw = text.encode("utf-8")
        buffer = C.create_string_buffer(raw)
        obj = self.send(self.cls("NSString"), "alloc")
        return self.send(obj, "initWithBytes:length:encoding:", C.c_void_p,
                         (C.c_void_p, C.c_size_t, C.c_size_t),
                         (C.cast(buffer, C.c_void_p), len(raw), 4))  # NSUTF8StringEncoding

    def text(self, obj):
        if not obj:
            return None
        description = self.send(obj, "description")
        raw = self.send(description, "UTF8String", C.c_char_p)
        require(raw is not None, "unreadable Objective-C description")
        return raw.decode("utf-8", errors="strict")

    def error(self, obj, seen=None):
        if not obj:
            return None
        seen = set() if seen is None else seen
        if obj in seen:
            return {"reference": hex(obj)}
        seen.add(obj)
        require(len(seen) <= 64, "NSError underlying graph exceeds diagnostic bound")
        info = self.send(obj, "userInfo")
        keys = self.send(info, "allKeys") if info else None
        values = {}
        underlying = []
        for index in range(self.integer(keys, "count") if keys else 0):
            key = self.send(keys, "objectAtIndex:", C.c_void_p, (C.c_size_t,), (index,))
            value = self.send(info, "objectForKey:", C.c_void_p, (C.c_void_p,), (key,))
            label = self.text(key)
            values[label] = self.text(value)
            if value and self.send(value, "isKindOfClass:", C.c_bool,
                                   (C.c_void_p,), (self.cls("NSError"),)):
                underlying.append({"key": label, "error": self.error(value, seen)})
            elif label in {"NSMultipleUnderlyingErrors", "NSMultipleUnderlyingErrorsKey"} and value:
                require(self.send(value, "isKindOfClass:", C.c_bool,
                                  (C.c_void_p,), (self.cls("NSArray"),)), "invalid underlying-error array")
                for number in range(self.integer(value, "count")):
                    item = self.send(value, "objectAtIndex:", C.c_void_p, (C.c_size_t,), (number,))
                    require(self.send(item, "isKindOfClass:", C.c_bool,
                                      (C.c_void_p,), (self.cls("NSError"),)), "invalid underlying-error object")
                    underlying.append({"key": label, "index": number, "error": self.error(item, seen)})
        return {"reference": hex(obj), "domain": self.text(self.send(obj, "domain")),
                "code": self.send(obj, "code", C.c_ssize_t),
                "localized_description": self.text(self.send(obj, "localizedDescription")),
                "localized_failure_reason": self.text(self.send(obj, "localizedFailureReason")),
                "description": self.text(obj), "user_info_description": self.text(info),
                "user_info": values, "underlying_errors": underlying}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--msl", type=Path, required=True)
    parser.add_argument("--expected-device-name", default="Apple Paravirtual device")
    parser.add_argument("--expected-registry-id", type=lambda value: int(value, 0))
    parser.add_argument("--optimization-level", choices=["default", "size"], default="default")
    parser.add_argument("--without-forced-inline", type=Path,
                        help="sealed annotation-only source control")
    args = parser.parse_args()
    require(sys.platform == "darwin" and platform.machine() == "arm64", "native arm64 Darwin required")
    payload = args.msl.read_bytes()
    require(len(payload) == MSL_BYTES and hashlib.sha256(payload).hexdigest() == MSL_SHA,
            "MSL differs from exact captured Initialize source")
    compiled = payload
    if args.without_forced_inline:
        compiled = args.without_forced_inline.read_bytes()
        require(compiled == unforced_msl(payload), "derived MSL changes more than forced-inline annotations")
    emit("input", path=str(args.msl.resolve()), bytes=len(payload), sha256=MSL_SHA,
         candidate_source=SOURCE, qualification_claimed=False,
         requested_optimization_level=args.optimization_level,
         compiled_path=str((args.without_forced_inline or args.msl).resolve()),
         compiled_bytes=len(compiled), compiled_sha256=hashlib.sha256(compiled).hexdigest(),
         removed_forced_inline_declarations=FORCED_INLINE_COUNT if args.without_forced_inline else 0,
         captured_source_modified=bool(args.without_forced_inline),
         route="synchronous native Metal error-returning methods; no command queue or dispatch")
    api = ObjC()
    pool = api.send(api.cls("NSAutoreleasePool"), "new")
    try:
        api.metal.MTLCopyAllDevices.argtypes = []
        api.metal.MTLCopyAllDevices.restype = C.c_void_p
        devices = api.metal.MTLCopyAllDevices()
        require(devices and api.integer(devices, "count") == 1, "expected exactly one Metal device")
        device = api.send(devices, "objectAtIndex:", C.c_void_p, (C.c_size_t,), (0,))
        name = api.text(api.send(device, "name"))
        registry = api.send(device, "registryID", C.c_uint64)
        require(name == args.expected_device_name and registry > 0, "Metal device name/registry identity differs")
        require(args.expected_registry_id is None or registry == args.expected_registry_id,
                "Metal registryID differs from explicit expected ID")
        second = api.metal.MTLCopyAllDevices()
        require(second and api.integer(second, "count") == 1, "Metal device inventory changed")
        selected = api.send(second, "objectAtIndex:", C.c_void_p, (C.c_size_t,), (0,))
        require(api.send(selected, "registryID", C.c_uint64) == registry
                and api.text(api.send(selected, "name")) == name, "Metal registry/name selection changed")
        device = selected
        emit("device", count=1, name=name, registry_id=registry,
             expected_registry_id=args.expected_registry_id,
             registry_scope="observed and re-enumerated this run; prior Vulkan bundle contains no Metal registryID",
             low_power=api.boolean(device, "isLowPower"),
             unified_memory=api.boolean(device, "hasUnifiedMemory"))
        options = api.send(api.cls("MTLCompileOptions"), "new")
        for selector, kind, value in [("setLanguageVersion:", C.c_size_t, (3 << 16) | 2),
                                      ("setFastMathEnabled:", C.c_bool, True),
                                      ("setPreserveInvariance:", C.c_bool, False)]:
            api.send(options, selector, None, (kind,), (value,))
        optimization = None
        requested_optimization = 1 if args.optimization_level == "size" else 0
        require(requested_optimization == 0 or api.supports(options, "setOptimizationLevel:"),
                "Metal library size optimization is unavailable")
        if api.supports(options, "setOptimizationLevel:"):
            api.send(options, "setOptimizationLevel:", None, (C.c_size_t,), (requested_optimization,))
            optimization = api.integer(options, "optimizationLevel")
            require(optimization == requested_optimization, "library optimization readback differs")
        require(api.integer(options, "languageVersion") == (3 << 16) | 2
                and api.boolean(options, "fastMathEnabled")
                and not api.boolean(options, "preserveInvariance"), "compile option readback differs")
        require(api.empty_collection(options, "preprocessorMacros") and api.empty_collection(options, "libraries"),
                "unexpected compile macros or dynamic libraries")
        emit("compile_options", language_version=api.integer(options, "languageVersion"),
             fast_math_enabled=api.boolean(options, "fastMathEnabled"),
             math_mode=api.integer(options, "mathMode") if api.supports(options, "mathMode") else None,
             math_fp_functions=api.integer(options, "mathFloatingPointFunctions") if api.supports(options, "mathFloatingPointFunctions") else None,
             preserve_invariance=False, optimization_level=optimization,
             changed_parameter={"forced_inline_declarations": {"baseline": FORCED_INLINE_COUNT, "diagnostic": 0}}
                 if args.without_forced_inline else
                 ({"optimization_level": {"baseline": 0, "diagnostic": 1}} if requested_optimization == 1 else None),
             comparison_control=("Exact captured MSL with the same " + args.optimization_level + " options")
                 if args.without_forced_inline else "Exact captured MSL with Default options",
             provenance="Original options with explicit diagnostic controls. Unforced source changes only annotations relative to the captured source with the same optimization level; no provider or numerical equivalence is claimed.")
        error = C.c_void_p()
        emit("phase_start", phase="library")
        library = api.send(device, "newLibraryWithSource:options:error:", C.c_void_p,
                           (C.c_void_p, C.c_void_p, C.POINTER(C.c_void_p)),
                           (api.string(compiled.decode("utf-8")), options, C.byref(error)))
        emit("phase_end", phase="library", success=bool(library), error=api.error(error.value))
        if not library:
            emit("terminal", pipeline_created=False, failed_stage="library",
                 qualification_claimed=False, execution_claimed=False)
            return 1
        emit("phase_start", phase="function")
        function = api.send(library, "newFunctionWithName:", C.c_void_p,
                            (C.c_void_p,), (api.string("main0"),))
        emit("phase_end", phase="function", success=bool(function))
        require(function, "main0 function missing")
        constants = api.send(function, "functionConstantsDictionary")
        require(not constants or api.integer(constants, "count") == 0, "unexpected function constants")
        descriptor = api.send(api.cls("MTLComputePipelineDescriptor"), "new")
        api.send(descriptor, "setComputeFunction:", None, (C.c_void_p,), (function,))
        api.send(descriptor, "setMaxTotalThreadsPerThreadgroup:", None, (C.c_size_t,), (1,))
        api.send(descriptor, "setThreadGroupSizeIsMultipleOfThreadExecutionWidth:", None, (C.c_bool,), (False,))
        # MoltenVK leaves these framework defaults untouched. Their getters may
        # return empty descriptor objects; nil identity is not their contract.
        observed = {
            "max_total_threads": api.integer(descriptor, "maxTotalThreadsPerThreadgroup"),
            "full_subgroups": api.boolean(descriptor, "threadGroupSizeIsMultipleOfThreadExecutionWidth"),
            "selected_function_matches": api.send(descriptor, "computeFunction") == function,
            "linked_functions": api.text(api.send(descriptor, "linkedFunctions")),
            "stage_input": api.text(api.send(descriptor, "stageInputDescriptor")),
            "binary_archives": api.text(api.send(descriptor, "binaryArchives")),
        }
        emit("pipeline_descriptor", entry_point="main0", local_size=[1, 1, 1],
             observed=observed, default_properties="untouched newly constructed framework defaults, as in MoltenVK",
             pipeline_options=0, description=api.text(descriptor))
        require(observed["max_total_threads"] == 1 and not observed["full_subgroups"]
                and observed["selected_function_matches"], "assigned compute descriptor fields differ")
        error = C.c_void_p()
        emit("phase_start", phase="pipeline")
        pipeline = api.send(device, "newComputePipelineStateWithDescriptor:options:reflection:error:", C.c_void_p,
                            (C.c_void_p, C.c_size_t, C.c_void_p, C.POINTER(C.c_void_p)),
                            (descriptor, 0, None, C.byref(error)))
        emit("phase_end", phase="pipeline", success=bool(pipeline), error=api.error(error.value))
        emit("terminal", pipeline_created=bool(pipeline), qualification_claimed=False,
             execution_claimed=False)
        return 0 if pipeline else 1
    finally:
        api.send(pool, "drain", None)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        emit("incomplete", error=f"{type(error).__name__}: {error}")
        raise SystemExit(1)
