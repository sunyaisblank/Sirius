"""Same original ABBA, counts and gates; both arms consume one genuine new host."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/benefit/prepare.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'ba3c35750f9e1b0b8e876b084d1f3f738ab34d1080709080dd5d4dc76915a440'
source = basis.read_text()
changes = [("doc(current / 'provider/linked-native-proof.json')",
  "doc(root / 'out/native-cpu-value-fact-stores/provider/linked-native-proof.json')"),
 ("build = doc(current / 'provider/provider-build-result.json')",
  "build = doc(root / 'out/native-cpu-value-fact-stores/provider/provider-build-result.json')"),
 ("assert pin(executable)['sha256'] == "
  "'7f2d1fa88953b954ab116233458a4cb962470f863fda9bbb3a1cad7d4a7d6a26'",
  "host = doc(root / 'out/native-cpu-owned-lower-increments/host-build-result.json')\n"
  "assert host['passed'] and host['new_cpp_compilations'] == 1\n"
  "assert host['all3modules_once_in_both_executables']\n"
  "owner = doc(root / 'out/native-cpu-owned-lower-increments/build/host-build/owner.json')\n"
  "assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']\n"
  "assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']\n"
  "assert str(executable) == str(root / host['executable'])\n"
  "assert pin(executable) == host['executable_identity']"),
 ('current B', 'current host/M'),
 ('same exact7f2d1fa8 host', 'same genuinely rebuilt host'),
 ('native B candidate', 'inline native R candidate')]
for before, after in changes:
    assert source.count(before) == 1, before
    source = source.replace(before, after)
source = source.replace('current host/M', 'current H host/AM').replace('Same genuinely rebuilt host/Mprovider', 'Same H host/AMprovider').replace('unchanged M property implementation', 'unchanged R property implementation').replace('exact M provider reused.', 'exact R provider used.')
source = source.replace('repaired inline native R candidate','AM complete-fact identity-store candidate').replace('current H host/AM','current H host/AM')
exec(compile(source, str(basis), "exec"))
