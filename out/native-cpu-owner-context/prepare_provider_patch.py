"""Reuse sealed prior authority; task-local paths/source-input binding only."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/prepare_provider_patch.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'e1e5c9757dd1f4d8272fd6b38b591b39fe0e8a7eb799cac0a0bc5a81514bb3f3'
source = basis.read_text()
changes = [("basis = root / 'out/native-cpu-slang-bridge'", "basis = root / 'out/native-cpu-active-prefix'"), ("root / 'bin/linux-gcc/native-cpu-active-prefix/typed-batch'", "root / 'bin/linux-gcc/native-cpu-owner-context/typed-batch'"), ("root / 'bin/linux-gcc/native-cpu-slang-bridge/typed-batch'", "root / 'bin/linux-gcc/native-cpu-active-prefix/typed-batch'"), ("('SIR-CPU-25.2.8A', 'SIR-CPU-25.2.8B')", "('SIR-CPU-25.2.8B', 'SIR-CPU-25.2.8C')"), ("('Sirius-native-cpu-slang-bridge-v1', 'Sirius-native-cpu-active-prefix-v2')", "('Sirius-native-cpu-active-prefix-v2', 'Sirius-native-cpu-owner-context-v3')"), ('Substantive activeprefix bridge repair; original shaderCPP/HPP/SPV/programs and host unchanged.', 'Privateownercontextlifetime repair only; originalCPP/SPV/programs/host and sharedBactiveprefix guard unchanged. Three derivednativeHPPs move only context declaration/globalParams assignment inside existinglanezero branch, retaining all64Groupcalls and exact ownerbody.'), ("'source_inputs': {'out/native-cpu-active-prefix/native_packet.hpp': pin((work / 'native_packet.hpp').read_bytes())}", "'source_inputs': {'out/native-cpu-active-prefix/native_packet.hpp': pin((root / 'out/native-cpu-active-prefix/native_packet.hpp').read_bytes()), 'out/native-cpu-owner-context/header-preparation.json': pin((work / 'header-preparation.json').read_bytes())}")]
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
