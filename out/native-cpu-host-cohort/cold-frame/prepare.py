"""Original frame settings/resources, genuine rebuilt host and unchanged M provider."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-simd-medium-tails/cold-frame/prepare.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'c2df3b7d5d1b04498ce3a0e2f50187f13a34d6d0a01e8c02c2d54b58bea4ba7e'
source = basis.read_text()
changes = [("EXE = R / 'bin/linux-gcc/native-cpu-dispatch-admission/frame/sirius_render_tests'",
  "EXE = R / 'bin/linux-gcc/native-cpu-host-preparation/frame/sirius_render_tests'"),
 ("ARCHIVE = R / 'bin/linux-gcc/native-cpu-row-ownership/matched/libsirius_backend.a'",
  "ARCHIVE = R / 'bin/linux-gcc/native-cpu-host-preparation/matched/libsirius_backend.a'"),
 ("assert doc(P / 'provider/linked-native-proof.json')['passed']",
  "assert doc(R / 'out/native-cpu-nir-dead-progress/provider/linked-native-proof.json')['passed']"),
 ("assert pin(EXE) == build['executable_identity'] == {\n"
  "    'bytes': 12469552, 'sha256': "
  "'a767656e1ba8faf8eb08b7ace82f2f076f8c7157d9f336ea2529137daecda1a1'}",
  "host = doc(R / 'out/native-cpu-host-preparation/host-build-result.json')\n"
  "assert host['passed'] and host['new_cpp_compilations'] == 1\n"
  "assert host['all3modules_once_in_both_executables']\n"
  "assert str(EXE) == str(R / host['frame_executable'])\n"
  "assert pin(EXE) == host['frame_executable_identity']"),
 ("assert pin(ARCHIVE) == recipe['backend'] == {\n"
  "    'bytes': 6871532, 'sha256': "
  "'07693aba054a519f1127220bcb3afd94e42d10f416c7c6358fda336af046f215'}",
  "assert str(ARCHIVE) == str(R / host['backend_archive'])\n"
  "assert pin(ARCHIVE) == host['backend_archive_identity']"),
 ("assert pin(HOST)['sha256'] == "
  "'7f2d1fa88953b954ab116233458a4cb962470f863fda9bbb3a1cad7d4a7d6a26'",
  "assert str(HOST) == str(R / host['executable'])\n"
  "assert pin(HOST) == host['executable_identity']"),
 ("assert recipe['argv'].count(str(ARCHIVE)) == 1\n"
  "assert recipe['argv'][recipe['argv'].index('-o') + 1] == str(EXE)",
  "assert recipe['argv'].count(str(R / "
  "'bin/linux-gcc/native-cpu-row-ownership/matched/libsirius_backend.a')) == 1\n"
  "assert recipe['argv'][recipe['argv'].index('-o') + 1] == str(R / "
  "'bin/linux-gcc/native-cpu-dispatch-admission/frame/sirius_render_tests')"),
 ("cache, destination = R / recovery['source'], R / recovery['destination']",
  "cache, destination = R / recovery['source'], EXE.parent / 'resources'"),
 ("write(W / 'reused-frame.json',", "write(W / 'frame-inputs.json',"),
 ("'historical_build': pin(F / 'build-result.json'), 'historical_link_recipe': pin(F / "
  "'link-recipe.json'),",
  "'historical_build': pin(F / 'build-result.json'), 'historical_link_recipe': pin(F / "
  "'link-recipe.json'),\n"
  "        'host_build': pin(R / 'out/native-cpu-host-preparation/host-build-result.json'),"),
 ("'backend_archive': pin(ARCHIVE), 'matched_7f2_host': pin(HOST), 'module_joins': joins,",
  "'backend_archive': pin(ARCHIVE), 'matched_host': pin(HOST), 'module_joins': joins,"),
 ("'resource_restore': pin(W / 'resource-restore.json'), 'new_links': 0, 'new_CPP_compilations': 0",
  "'resource_restore': pin(W / 'resource-restore.json'), 'new_links': 2, 'new_CPP_compilations': "
  '1'),
 ("reuse = doc(W / 'reused-frame.json')", "reuse = doc(W / 'frame-inputs.json')"),
 ("assert reuse['passed'] and reuse['new_links'] == reuse['new_CPP_compilations'] == 0",
  "assert reuse['passed'] and reuse['new_links'] == 2 and reuse['new_CPP_compilations'] == 1\n"
  "    assert reuse['host_build'] == pin(R / "
  "'out/native-cpu-host-preparation/host-build-result.json')"),
 ("'control-8/native-observation.json', 'control-8/inputs-before.json', "
  "'provider/linked-native-proof.json'",
  "'control-8/native-observation.json', 'control-8/inputs-before.json'"),
 ('inputs |= {path for path in W.iterdir() if path.is_file()}',
  'inputs |= {path for path in W.iterdir() if path.is_file()}\n'
  "    inputs |= {R / 'out/native-cpu-host-preparation/host-build-result.json',\n"
  "               R / 'out/native-cpu-nir-dead-progress/provider/linked-native-proof.json',\n"
  "               R / 'out/native-cpu-runtime-observe/all_threads.py'}"),
 ('renders/native-cpu-simd-medium-tails-', 'renders/native-cpu-host-cohort-'),
 ("assert action['argv'] == ['/usr/bin/cmake', '-DSIRIUS_RENDER_TESTS=' + str(EXE), '-P',\n"
  "                              str(R / 'tests/operational/vulkan_runtime_profile.cmake')]",
  "assert action['argv'] == ['/usr/bin/cmake', '-DSIRIUS_RENDER_TESTS=' + str(R / "
  "'bin/linux-gcc/native-cpu-dispatch-admission/frame/sirius_render_tests'), '-P',\n"
  "                              str(R / 'tests/operational/vulkan_runtime_profile.cmake')]\n"
  "    action['argv'][1] = '-DSIRIUS_RENDER_TESTS=' + str(EXE)"),
 ('Exact preserved F frameELF and matched7f2 backend reused; no fresh link or C++ compilation.',
  'Genuine new host/frame ELF with one changed C++ unit and two audited links; exact M provider '
  'reused.'),
 ("print(json.dumps({'restored_names': 32, 'resource_payload_bytes_copied': 0, 'new_links': 0, "
  "'executed': False}))",
  "print(json.dumps({'restored_names': 32, 'resource_payload_bytes_copied': 0, "
  "'host_links_already_completed': 2, 'executed': False}))"),
 ("print(json.dumps({'prepared_original_frame': True, 'new_links': 0, 'executed': False}))",
  "print(json.dumps({'prepared_original_frame': True, 'host_links_already_completed': 2, "
  "'executed': False}))")]
for before, after in changes:
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), "exec"))
