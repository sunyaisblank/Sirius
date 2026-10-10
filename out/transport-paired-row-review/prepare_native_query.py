"""Prepare one bounded ordinary AMD compiler query, with zero dispatches."""
from pathlib import Path, PureWindowsPath
import hashlib
import json
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
TASK = Path(__file__).resolve().parent
OLD = ROOT / 'attestations/native-vulkan/3045067/ordinary-pipeline-generated-code'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def win(path):
    return r'\\wsl.localhost\Ubuntu' + str(path.resolve()).replace('/', '\\')

revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
assert revision == '6845c8150ebf0df6992fe2609a4e9323c7f0612b'
assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT)
compiled = json.loads((TASK / 'compile-result.json').read_text())
assert compiled['compile_success'] and compiled['all_selected_inputs_unchanged']
module = TASK / 'paired-transport.spv'
assert sha(module) == '721be7625f8eebbf722adf3e79a89050bffe9bf9b7260b8f5120b0a85ccd7f27'
query = TASK / 'native-query'
query.mkdir(exist_ok=False)
provider = json.loads((OLD / 'run-query/driver-modules.json').read_text(encoding='utf-8-sig'))
provider_bindings = []
for item in provider:
    path = PureWindowsPath(item['path'])
    host = Path('/mnt/c').joinpath(*path.parts[1:])
    assert host.stat().st_size == item['bytes'] and sha(host) == item['sha256']
    provider_bindings.append({'path': str(host), 'bytes': host.stat().st_size, 'sha256': sha(host)})
source = (OLD / 'native_query.cs').read_text()
replacements = {
    'string[] paths = new string[] { transportPath, endpointPath };':
    'string[] paths = new string[] { transportPath };',
    'string[] stages = new string[] { "Transport", "Endpoint" };':
    'string[] stages = new string[] { "TransportPaired" };',
    'string[] expectedHashes = new string[] { "b52e6373f6b30aac6c7fe877f6e07f3c47c7ac9eafa1671840eacd19960d7cb1", "785d94c53839693af8a0c752b8d7f0f1bd44789df17f3ab4eef79d3684faa559" };':
    'string[] expectedHashes = new string[] { "' + sha(module) + '" };',
    '"pipelines_created", 2, "dispatches", 0': '"pipelines_created", 1, "dispatches", 0',
    'Unchanged production module/features/layout/flags with AMD tooling extension added. Static generated-code observations; no dispatch, timing comparison, numerical or full qualification verdict.':
    'Isolated unadopted paired-row Transport module; production features/layout/flags, tooling extension added. Static generated-code only; zero dispatches, no numerical, speed, or qualification verdict.',
}
for before, after in replacements.items():
    assert source.count(before) == 1, before
    source = source.replace(before, after)
(query / 'native_query.cs').write_text(source)
child = (OLD / 'query-child.ps1').read_text()
child = child.replace(win(ROOT / 'out/native-ordinary-shader-info'), win(query))
child = child.replace(win(query / 'query-payloads/transport.spv'), win(module))
child = child.replace(win(query / 'query-payloads/endpoint.spv'), win(module))
child = child.replace(win(query / 'query-payloads'), win(query / 'payloads'))
(query / 'query-child.ps1').write_text(child)
owner = (OLD / 'query-owner.ps1').read_text()
owner = owner.replace('$result.pipelines_created -ne 2 -or @($result.observations).Count -ne 2',
                      '$result.pipelines_created -ne 1 -or @($result.observations).Count -ne 1')
owner = owner.replace('3045067d4c29da744059798e175b97942a319f84', revision)
(query / 'query-owner.ps1').write_text(owner)
for name in ('query-native-temp', 'payloads'):
    (query / name).mkdir()
source_paths = [OLD / name for name in ('native_query.cs', 'query-child.ps1', 'query-owner.ps1')]
source_paths += [ROOT / 'src/sirius/backend/vulkan/vulkan_device.cpp', ROOT / 'src/sirius/backend/retained_compute.cpp',
                 ROOT / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h',
                 TASK / 'prototype.py', TASK / 'prepare_native_query.py', module,
                 TASK / 'kernels/retained_transport.slang', TASK / 'kernels/retained_scratch.slang']
source_paths += [query / name for name in ('native_query.cs', 'query-child.ps1', 'query-owner.ps1')]
source_paths += [Path('/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe')]
bindings = [{'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha(path)} for path in source_paths]
record = {'revision': revision, 'prototype_only': True, 'device_dispatches_authorized': 0,
          'query_pipelines': 1, 'pipeline_flags': 0, 'stage_flags': 0, 'empty_application_cache': True,
          'driver_cold_claimed': False, 'seconds_guard': 90, 'rss_bytes_guard': 4294967296,
          'provider_bindings': provider_bindings, 'input_bindings': bindings,
          'scope': 'Bounded static AMD compiler/query only. Original scientific/render governors unchanged; no device dispatch or numerical/performance acceptance.'}
(query / 'preparation.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'query': str(query), 'selected_inputs': len(bindings), 'provider_inputs': len(provider_bindings), 'module_sha256': sha(module)}, indent=2))
