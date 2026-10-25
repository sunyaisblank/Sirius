"""Bind original C8 to real P host/new R provider and fresh Transport coverage."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/prepare.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '6a414a90ea995c5cfb93d3c132c91397fe9203160d930f80125d2f8785ed85cb'
source = basis.read_text()
changes = [("native = root / 'out/native-cpu-active-prefix'", "native = root / 'out/native-cpu-reachable-tag-writers-cohort'"),
 ("build = doc(native / 'provider/provider-build-result.json')",
  "build = doc(root / 'out/native-cpu-reachable-tag-writers/provider/provider-build-result.json')"),
 ("compiled = doc(native / 'batch-compile-result.json')",
  "compiled = doc(root / 'out/native-cpu-scalar-zero/cpp-compile-result.json')"),
 ("doc(native / 'provider/linked-native-proof.json')",
  "doc(root / 'out/native-cpu-reachable-tag-writers/provider/linked-native-proof.json')"),
 ("executable = root / 'bin/linux-gcc/native-cpu-row-ownership/matched/audited_backend_tests'",
  "executable = root / 'bin/linux-gcc/native-cpu-owned-lower-increments/matched/audited_backend_tests'\n"
  "host = doc(root / 'out/native-cpu-owned-lower-increments/host-build-result.json')\n"
  "assert host['passed'] and host['new_cpp_compilations'] == 1\n"
  "assert host['all3modules_once_in_both_executables']\n"
  "owner = doc(root / 'out/native-cpu-owned-lower-increments/build/host-build/owner.json')\n"
  "assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']\n"
  "assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']\n"
  "assert str(executable) == str(root / host['executable'])\n"
  "assert pin(executable) == host['executable_identity']"),
 ("expected.update(compiled['source_inputs'])",
  "expected.update(compiled['source_inputs'])\n"
  "expected.update(host['source_inputs'])\n"
  'expected.update(doc(root / '
  "'out/native-cpu-reachable-tag-writers/provider/source-contract.json')['source_inputs'])\n"
  "for name in ('host-build-result.json', 'build/host-build/owner.json'):\n"
  "    path = root / 'out/native-cpu-owned-lower-increments' / name\n"
  '    expected[str(path)] = pin(path)\n'
  "expected[str(executable)] = host['executable_identity']\n"
  "expected[str(root / host['backend_archive'])] = host['backend_archive_identity']\n"
  "expected[str(root / host['frame_executable'])] = host['frame_executable_identity']"),
 ("assert pin(executable)['sha256'] == "
  "'7f2d1fa88953b954ab116233458a4cb962470f863fda9bbb3a1cad7d4a7d6a26'",
  "assert pin(executable) == host['executable_identity']"),
 ("(native / 'provider').iterdir()",
  "(root / 'out/native-cpu-reachable-tag-writers/provider').iterdir()"),
 ('unchanged complete-row matched ELF7f2d1fa8 and nativeCPU active-prefix providerB',
  'genuinely rebuilt complete-row host with exact normal spacing/redundant active zero-store '
  'removal and nativeCPU providerAL with Transport worker occupancy scheduling'),
 ('First actual native packet/return perstage observed passively;',
  'First actual rows/pure/SIMD8 entries and returns perstage observed passively;')]
for before, after in changes:
    assert source.count(before) == 1, before
    source = source.replace(before, after)
before = "files = set(expected)"
after = """compatibility_namespace = {}
compatibility_source = work / 'compatibility.py'
exec(compile(compatibility_source.read_text(), str(compatibility_source), 'exec'), compatibility_namespace)
component = compatibility_namespace['check'](root)
expected.update(component['inputs'])
files = set(expected)"""
assert source.count(before) == 1
source = source.replace(before, after)
source = source.replace('nativeCPU providerAL with Transport worker occupancy scheduling',
    'nativeCPU providerAL with unchanged Q Transport scheduling and identical U scalar3, AL raw3 with corrected reachable tag invalidation/refresh, and exact AK arithmetic/adapter including its six empty-work Divide/Sqrt returns, final native NIR optimization deferred with full fallback retained; fresh Transport compatibility required')
source = source.replace('First actual rows/pure/SIMD8 entries and returns perstage observed passively;',
    'First actual rows/pure returns and policy-selected implementations observed passively; changed Transport SIMD has fresh full38-case compatibility;')
before = 'for name, value in expected.items():'
assert source.count(before) == 1
source = source.replace(before, """runtime_namespace = {}
runtime_source = work / 'runtime_transition.py'
exec(compile(runtime_source.read_text(), str(runtime_source), 'exec'), runtime_namespace)
runtime_namespace['apply'](root, work, expected)
for name, value in expected.items():""")
exec(compile(source, str(basis), "exec"))
