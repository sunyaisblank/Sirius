"""Reuse reviewed bounded native lifecycle; change only this batch's bindings."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
OLD=ROOT/'attestations/native-vulkan/24a8e14/native-rhs-layer-clock-execution/files'
review=json.loads((OLD/'controller-review.json').read_text());assert review['accepted']
def seal(p):
    data=p.read_bytes();return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
records={}
for name in ('run_native_control.py','run_shader_clock_control.ps1','emergency_native_cleanup.ps1','run_native_retirement.py','verify_native_retirement.ps1','prepare_native_stage.py'):
    p=OLD/name;assert seal(p)==review['controller_seals'][name]
    text=p.read_text()
    if name=='run_native_control.py':
        text=text.replace('Original numerical control with shader-clock sidecars', 'Original numerical control with word-value sidecars')
        text=text.replace('Root-owned original numerical control with shader-clock sidecars', 'Root-owned original numerical control with word-value sidecars')
    if name=='run_shader_clock_control.ps1':
        text=text.replace('Uncalibrated diagnostic clock sidecars cleared and captured around original synchronous commands.', 'Diagnostic word-value sidecars cleared and complete original input/output snapshots captured after original synchronous commands.')
    if name=='prepare_native_stage.py':
        text=text.replace("assert HEAD == '24a8e145d1da5dca75a4fa26ec8e639020eecf37'", "SOURCE = json.loads((WORK / 'source.json').read_text())\nassert HEAD == SOURCE['source_revision']")
        text=text.replace("STAGE = ROOT / 'bin/windows-msvc/native-24a8e14-rhs-layer-clock'", "STAGE = ROOT / SOURCE['stage_root']")
        old="path = (WORK / 'phase/retained_transport.spv' if name == 'retained_transport.spv' else ROOT / 'attestations/native-vulkan/2053ef7/native-shader-clock-admission/files/phase/retained_endpoint.spv')"
        assert text.count(old)==1; text=text.replace(old,"path = WORK / 'kernels' / name\n    stage = 'transport' if name == 'retained_transport.spv' else 'endpoint'\n    reviewed = json.loads((WORK / 'compiled-review.json').read_text())\n    assert reviewed['accepted'] and seal(path) == reviewed['modules'][stage]")
        start=text.index("for path in ('src/sirius/backend/vulkan/vulkan_device.cpp'")
        end=text.index('plan = {',start)
        text=text[:start]+'''tracked = subprocess.check_output(['git', 'ls-files', '-z']).split(b'\\0')
for raw in tracked:
    if raw:
        path = raw.decode(); inputs[path] = seal(ROOT / path)
for path in list(WORK.iterdir()) + list((WORK / 'kernels').iterdir()):
    if path.is_file() and path.name not in ('execution-plan.json', 'controller-review.json'):
        inputs[str(path.relative_to(ROOT))] = seal(path)
assert json.loads((WORK / 'host-source-review.json').read_text())['accepted']
assert json.loads((WORK / 'ci.json').read_text())['accepted']
assert json.loads((WORK / 'artifact-extraction.json').read_text())['accepted']
''' + text[end:]
        text=text.replace('uncalibrated native-default shader-clock phase/layer sidecars and exact input snapshots', 'native-default word-value sidecars and complete actual input/output snapshots')
        text=text.replace("len(inputs) - 9", "sum(p.is_file() for p in STAGE.rglob('*'))")
    target=WORK/name;assert not target.exists();target.write_text(text)
    records[name]={'basis':str(p.relative_to(ROOT)),'basis_seal':seal(p),'current_seal':seal(target),'byte_identical':p.read_bytes()==target.read_bytes()}
(WORK/'native-source-preparation.json').write_text(json.dumps({'records':records,'scope':'Unexecuted private native lifecycle preparation. Old320s4GiB/400s/40s controls and8MiB original allocation retained. Fresh source/CI/export/module/provider checks required.'},indent=2)+'\n')
print(json.dumps({k:v['current_seal'] for k,v in records.items()},indent=2))
