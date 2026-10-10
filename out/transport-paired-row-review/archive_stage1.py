"""Preserve the completed isolated prerequisite without qualification claims."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
TASK = Path(__file__).resolve().parent
DESTINATION = ROOT / 'attestations/native-vulkan/6845c81/paired-transport-source-numerical'
EXCLUDED = {'numerical-prepare-initial', 'numerical-prepare-reviewed-first',
            'native-temp', 'query-owner-temp', '__pycache__'}


def binding(path):
    data = path.read_bytes()
    return {'path': str(path), 'bytes': len(data),
            'sha256': hashlib.sha256(data).hexdigest()}


assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip() == '6845c8150ebf0df6992fe2609a4e9323c7f0612b'
assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT)
assert not DESTINATION.exists()
result = json.loads((TASK / 'numerical/reference-check.json').read_text())
assert result['pass'] and result['dispatches'] == 54 and result['late_refused_rows'] == 14
bridge = json.loads((TASK / 'numerical/bridge-query.json').read_text())
assert bridge['returncode'] == 0 and bridge['all_selected_inputs_unchanged']
for item in bridge['input_bindings']:
    assert binding(Path(item['path'])) == item, item['path']
audit = json.loads((TASK / 'numerical/terminal-query-audit.json').read_text(encoding='utf-8-sig'))
assert audit['owned_births_absent'] and audit['all_two_module_seals_match']
assert not audit['consumer_scan']['matches'] and not audit['consumer_scan']['accessible_compilers']
payloads = []


def copy(source, destination):
    before = binding(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    copied = binding(destination)
    assert (before['bytes'], before['sha256']) == (copied['bytes'], copied['sha256'])
    payloads.append({'path': str(destination.relative_to(DESTINATION)),
                     'bytes': copied['bytes'], 'sha256': copied['sha256']})


for source in sorted(TASK.rglob('*')):
    relative = source.relative_to(TASK)
    if not source.is_file() or EXCLUDED.intersection(relative.parts):
        continue
    if source.name in ('archive-protection.json', 'source-git.index'):
        continue
    copy(source, DESTINATION / 'diagnostic' / relative)
# Preserve selected canonical source and the actual pre-integration generated
# header. Runtime providers/tool executables remain pinned rather than copied.
for item in bridge['input_bindings']:
    source = Path(item['path'])
    if not source.is_relative_to(ROOT) or source.is_relative_to(TASK):
        continue
    if source.suffix not in ('.py', '.cpp', '.h', '.slang', '.json', '.ps1', '.cs'):
        continue
    copy(source, DESTINATION / 'source-inputs' / source.relative_to(ROOT))
for source in (ROOT / 'scripts/build-retained-kernels.py',
               ROOT / 'src/sirius/kernels/retained_program.py',
               ROOT / 'src/sirius/kernels/retained_transport.slang',
               ROOT / 'src/sirius/kernels/retained_scratch.slang',
               ROOT / 'src/sirius/backend/retained_compute.h'):
    target = DESTINATION / 'source-inputs' / source.relative_to(ROOT)
    if not target.exists():
        copy(source, target)
manifest = {'baseline_revision': bridge['revision'], 'payloads': payloads,
            'payload_count': len(payloads), 'payload_bytes': sum(i['bytes'] for i in payloads),
            'scope': 'Isolated source/static native resource and finite raw/reference prerequisite. Previous stopped launcher and reachability failure preserved. No host integration, timing gain, original frame, full scientific/platform/release or adoption verdict.'}
(DESTINATION / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
receipt = {'archive': str(DESTINATION), 'manifest': binding(DESTINATION / 'manifest.json'),
           'payload_count': len(payloads), 'payload_bytes': manifest['payload_bytes']}
(TASK / 'archive-protection.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt, indent=2))
