from pathlib import Path
import hashlib
import json
import re
import shutil
import subprocess

root = Path.cwd()
work = root / 'out/macos-provider-metadata-control'
binary = root / 'bin/linux-gcc/macos-provider-metadata-control'
assert not work.exists() and not binary.exists()
work.mkdir()
binary.mkdir()
payload = root / 'out/macos-moltenvk-identity-fixed/payload'
device = 'sources/moltenvk-patched/MoltenVK/MoltenVK/GPUObjects/MVKDevice.mm'
content = json.loads((payload / 'content-id.json').read_text())
old = (payload / 'historical/original' / device).read_text()
new = (payload / device).read_text()
old_statement = re.findall(r'\tstrcpy\(supportedProps12\.driverInfo,.*?;', old, re.S)
new_statement = re.findall(r'\tstd::snprintf\(supportedProps12\.driverInfo,.*?;', new, re.S)
assert len(old_statement) == len(new_statement) == 1
header = '#include <mvk_private_api.h>\n#include <cstdio>\n#include <cstring>\n'
header += '#define MVK_DIAGNOSTIC_CONTENT_ID "' + content['content_id'] + '"\n'
tail = '\n  std::puts(supportedProps12.driverInfo);\n  return std::strcmp(supportedProps12.driverInfo, "' + content['driver_info_expected'] + '") != 0;\n}\n'
for name, statement in [('original', old_statement[0]), ('corrected', new_statement[0])]:
    (work / (name + '.cpp')).write_text(header + 'int main() {\n  VkPhysicalDeviceVulkan12Properties supportedProps12{};\n  static_assert(sizeof(supportedProps12.driverInfo) == 256);\n' + statement + tail)
for name in ('run_owned.py', 'owned_guard.py'):
    shutil.copy2(root / 'out/science-p4-cpu-remainder' / name, work / name)
def pin(path):
    raw = path.read_bytes()
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}
inputs = [payload / device, payload / 'historical/original' / device,
          payload / 'content-id.json', payload / 'input-freeze.json',
          payload / 'sources/moltenvk-patched/MoltenVK/MoltenVK/API/mvk_private_api.h',
          work / 'original.cpp', work / 'corrected.cpp', Path(__file__).resolve()]
inputs += [p for p in (payload / 'sources/Vulkan-Headers/include').rglob('*') if p.is_file()]
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
common = {'revision': revision, 'kind': 'check',
          'limits': {'seconds': 60, 'rss_kib': 4194304, 'reserve_kib': 2097152},
          'scope': 'Pinned actual-header diagnostic metadata expression only; Linux GCC syntax/output control, no native MoltenVK or numerical acceptance.',
          'frozen_inputs': [str(p.relative_to(root)) for p in inputs]}
flags = ['/usr/bin/g++-14', '-std=c++2c', '-Wall', '-Wextra', '-Werror',
         '-isystem', str(payload / 'sources/Vulkan-Headers/include'),
         '-isystem', str(payload / 'sources/moltenvk-patched/MoltenVK/MoltenVK/API')]
actions = {}
for name, expected in [('original', 1), ('corrected', 0)]:
    actions[name + '-compile'] = {**common, 'expected_exit': expected,
        'argv': [*flags, str(work / (name + '.cpp')), '-o', str(binary / name)]}
actions['corrected-run'] = {**common, 'expected_exit': 0, 'argv': [str(binary / 'corrected')]}
(work / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
(work / 'preparation.json').write_text(json.dumps({'runner_source': revision, 'content_id': content['content_id'],
    'expected_stdout': content['driver_info_expected'] + '\n', 'no_native_acceptance_claimed': True,
    'original_statement': old_statement[0], 'corrected_statement': new_statement[0],
    'inputs': {str(p.relative_to(root)): pin(p) for p in inputs}}, indent=2) + '\n')
print(json.dumps({'work': str(work), 'pins': len(inputs), 'content_id': content['content_id']}))
