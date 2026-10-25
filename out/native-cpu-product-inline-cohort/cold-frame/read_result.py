"""Original frame verdict and module/owner checks with genuine changed-host build proof."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-dispatch-admission/cold-frame/read_result.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '04b2ff93f78aa150dd4921cbe28aa458b21860d07692bb89a9f44d99b4a73dcf'
source = basis.read_text()
changes = [("current = root / 'out/native-cpu-dispatch-admission'",
  "current = root / 'out/native-cpu-product-inline-cohort'"),
 ("build = doc(work / 'build-result.json')\n"
  "assert build['passed'] and build['new_CPP_compilations'] == 0",
  "build = doc(root / 'out/native-cpu-upload-aware-staging/host-build-result.json')\n"
  "assert build['passed'] and build['new_cpp_compilations'] == 1\n"
  "assert build['all3modules_once_in_both_executables']\n"
  "restored = doc(work / 'resource-restore.json')\n"
  "frame_inputs = doc(work / 'frame-inputs.json')\n"
  "assert frame_inputs['passed'] and frame_inputs['new_links'] == 2 and "
  "frame_inputs['new_CPP_compilations'] == 1\n"
  "assert frame_inputs['host_build'] == pin(root / "
  "'out/native-cpu-upload-aware-staging/host-build-result.json')\n"
  "assert restored['status'] == 'complete' and restored['resource_payload_bytes_copied'] == 0\n"
  "assert frame_inputs['resource_restore'] == pin(work / 'resource-restore.json')\n"
  "cache, destination = root / restored['cache'], root / restored['destination']\n"
  "assert {str(path.relative_to(destination)) for path in destination.rglob('*') if "
  "path.is_file()} == set(restored['files'])\n"
  "for name, value in restored['files'].items():\n"
  '    assert (cache / name).samefile(destination / name)\n'
  '    assert pin(cache / name) == pin(destination / name) == value\n'),
 ("executable = root / build['executable']\nassert pin(executable) == build['executable_identity']",
  "executable = root / build['frame_executable']\n"
  "assert pin(executable) == build['frame_executable_identity']"),
 ("'build': pin(work / 'build-result.json'),",
  "'build': pin(root / 'out/native-cpu-upload-aware-staging/host-build-result.json'), 'frame_inputs': "
  "pin(work / 'frame-inputs.json'),"),
 ("for kind, row in build['module_joins'].items():\n    original = modules[kind]",
  "for kind, row in build['module_joins'].items():\n"
  "    assert row['executable_complete_occurrences'] == "
  "row['frame_executable_complete_occurrences'] == 1\n"
  "    row = {key: row[key] for key in ('module', 'bytes', 'sha256', 'complete_occurrences')}\n"
  '    original = modules[kind]')]
for before, after in changes:
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), "exec"))
