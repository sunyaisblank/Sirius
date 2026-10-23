"""Reuse the original complete-code observer with exact new provider/module paths."""
from pathlib import Path
import hashlib
P = Path(__file__).resolve().parents[2] / 'lavapipe-acqrel-ordering/ordinary/capture.py'
assert hashlib.sha256(P.read_bytes()).hexdigest() == '803758e9ce81b34330283c9e81a99bff92647ad176e60a10e1679086b82a40dd'
source = P.read_text()
source = source.replace("bin/linux-gcc/lavapipe-acqrel-ordering/src/gallium/targets/lavapipe/libvulkan_lvp.so", "bin/linux-gcc/lavapipe-intergroup-row-owner/provider/libvulkan_lvp.so")
source = source.replace("out/native-original-input-ownership/compiler-result.json", "out/native-cpu-row-ownership/compiler-result.json")
source = source.replace("Exact prior admitted complete module, new guarded acquire-release provider", "Exact previously audited row-owner module, new intergroup SIMD provider")
exec(compile(source, str(P), 'exec'))
