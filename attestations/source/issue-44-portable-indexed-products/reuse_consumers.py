"""Recheck exact original consumer and provider inputs without another baseline build."""
from pathlib import Path
import hashlib
import json
import resource
import shutil
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024 * 1024, 64 * 1024 * 1024))
root = Path.cwd()
work = Path(__file__).resolve().parent
build = root / 'bin/linux-gcc/portable-indexed-limb-accumulation'
archive = root / 'attestations/software-vulkan/d910f3c/portable-transport-limb-accumulation/files'
old = archive / 'out/portable-transport-limb-accumulation'
oldbuild = archive / 'bin/linux-gcc/portable-transport-limb-accumulation'
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
reuse = json.loads((old / 'baseline-reuse.json').read_text())
for item in reuse['actual_mutable_link_archives_tools_unchanged']:
    path = Path(item['path'])
    assert path.stat().st_size == item['bytes'] and sha(path) == item['sha256'], path
assert sha(oldbuild / 'original_baseline_controls') == reuse['baseline_executable_sha256']
assert sha(oldbuild / 'vulkan.o') == reuse['same_exact_vulkan_provider_object']
assert sha(build / 'vulkan-device-v2.o') == reuse['same_exact_vulkan_provider_object']
for name in ('original_baseline_controls', 'baseline.map', 'vulkan.o'):
    assert not (build / name).exists()
    shutil.copy2(oldbuild / name, build / name)
source = old / 'original-controls/candidate'
target = work / 'original-controls/candidate'
assert not target.exists()
shutil.copytree(source, target)
metadata = work / 'control-contracts'
assert not metadata.exists()
metadata.mkdir()
for name in ('software_contracts.py', 'expected-readback-layout.json'):
    shutil.copy2(old / 'control-contracts' / name, metadata / name)
reuse['preparation_current_revision'] = __import__('subprocess').check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
reuse['reuse_origin_archive'] = str(archive.parent.relative_to(root))
reuse['reuse_scope'] = 'Exact archived stock consumer/map and current Vulkan provider, with all recorded shared tools/archives rechecked. Fresh candidate objects/link and full 80 payload joins remain required; no current runtime or benefit verdict.'
(work / 'baseline-reuse.json').write_text(json.dumps(reuse, indent=2) + '\n')
print(json.dumps({'passed': True, 'shared_inputs_exact': len(reuse['actual_mutable_link_archives_tools_unchanged']),
                  'baseline_executable_exact': True, 'provider_object_exact': True}))
