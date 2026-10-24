"""Bind original C8 to one genuine changed host and unchanged finite-qualified M."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/prepare.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '6a414a90ea995c5cfb93d3c132c91397fe9203160d930f80125d2f8785ed85cb'
source = basis.read_text()
changes = [("native = root / 'out/native-cpu-active-prefix'", "native = root / 'out/native-cpu-host-cohort'"),
 ("build = doc(native / 'provider/provider-build-result.json')",
  "build = doc(root / 'out/native-cpu-nir-dead-progress/provider/provider-build-result.json')"),
 ("compiled = doc(native / 'batch-compile-result.json')",
  "compiled = doc(root / 'out/native-cpu-simd-endpoint-regions/ispc-uniform-compile-result.json')"),
 ("doc(native / 'provider/linked-native-proof.json')",
  "doc(root / 'out/native-cpu-nir-dead-progress/provider/linked-native-proof.json')"),
 ("executable = root / 'bin/linux-gcc/native-cpu-row-ownership/matched/audited_backend_tests'",
  "executable = root / 'bin/linux-gcc/native-cpu-host-preparation/matched/audited_backend_tests'\n"
  "host = doc(root / 'out/native-cpu-host-preparation/host-build-result.json')\n"
  "assert host['passed'] and host['new_cpp_compilations'] == 1\n"
  "assert host['all3modules_once_in_both_executables']\n"
  "owner = doc(root / 'out/native-cpu-host-preparation/build/host-build/owner.json')\n"
  "assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']\n"
  "assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']\n"
  "assert str(executable) == str(root / host['executable'])\n"
  "assert pin(executable) == host['executable_identity']"),
 ("expected.update(compiled['source_inputs'])",
  "expected.update(compiled['source_inputs'])\n"
  "expected.update(host['source_inputs'])\n"
  'expected.update(doc(root / '
  "'out/native-cpu-nir-dead-progress/provider/source-contract.json')['source_inputs'])\n"
  "for name in ('host-build-result.json', 'build/host-build/owner.json'):\n"
  "    path = root / 'out/native-cpu-host-preparation' / name\n"
  '    expected[str(path)] = pin(path)\n'
  "expected[str(executable)] = host['executable_identity']\n"
  "expected[str(root / host['backend_archive'])] = host['backend_archive_identity']\n"
  "expected[str(root / host['frame_executable'])] = host['frame_executable_identity']"),
 ("assert pin(executable)['sha256'] == "
  "'7f2d1fa88953b954ab116233458a4cb962470f863fda9bbb3a1cad7d4a7d6a26'",
  "assert pin(executable) == host['executable_identity']"),
 ("(native / 'provider').iterdir()",
  "(root / 'out/native-cpu-nir-dead-progress/provider').iterdir()"),
 ('unchanged complete-row matched ELF7f2d1fa8 and nativeCPU active-prefix providerB',
  'genuinely rebuilt complete-row host with exact normal spacing/redundant active zero-store '
  'removal and unchanged nativeCPU providerM'),
 ('First actual native packet/return perstage observed passively;',
  'First actual rows/pure/SIMD8 entries and returns perstage observed passively;')]
for before, after in changes:
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), "exec"))
