"""Reuse the original provider seam; emit only three reviewed file changes."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-active-prefix/prepare_provider_patch.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'e1e5c9757dd1f4d8272fd6b38b591b39fe0e8a7eb799cac0a0bc5a81514bb3f3'
source = basis.read_text()
changes = [
    ("basis = root / 'out/native-cpu-slang-bridge'", "basis = root / 'out/native-cpu-owner-context'"),
    ("root / 'bin/linux-gcc/native-cpu-active-prefix/typed-batch'", "root / 'bin/linux-gcc/native-cpu-static-program/typed-batch'"),
    ("root / 'bin/linux-gcc/native-cpu-slang-bridge/typed-batch'", "root / 'bin/linux-gcc/native-cpu-owner-context/typed-batch'"),
    ("('SIR-CPU-25.2.8A', 'SIR-CPU-25.2.8B')", "('SIR-CPU-25.2.8C', 'SIR-CPU-25.2.8D')"),
    ("('Sirius-native-cpu-slang-bridge-v1', 'Sirius-native-cpu-active-prefix-v2')", "('Sirius-native-cpu-owner-context-v3', 'Sirius-native-cpu-static-program-v4')"),
    ('Substantive activeprefix bridge repair; original shaderCPP/HPP/SPV/programs and host unchanged.',
     'Complete canonical native lowering: all9380instructions/399layers, original arithmetic and fallback/layer/refusal timing retained. OriginalCPP/SPV/programs/host/sharedBguard unchanged; Cownercontext retained.'),
    ('uniqueprivateBcache/UUID', 'uniqueprivateDcache/UUID'),
    ("'source_inputs': {'out/native-cpu-active-prefix/native_packet.hpp': pin((work / 'native_packet.hpp').read_bytes())}",
     "'source_inputs': {'out/native-cpu-active-prefix/native_packet.hpp': pin((root / 'out/native-cpu-active-prefix/native_packet.hpp').read_bytes()), 'out/native-cpu-static-program/header-preparation.json': pin((work / 'header-preparation.json').read_bytes())}"),
]
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
