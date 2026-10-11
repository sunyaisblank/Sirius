"""Reuse accepted native controllers and adapt only this batch's source bindings."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
OLD = ROOT / 'attestations/native-vulkan/8f1118d/native-shader-clock-execution/files'
HEAD = '24a8e145d1da5dca75a4fa26ec8e639020eecf37'
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == HEAD
assert subprocess.check_output(['git', 'rev-parse', '@{u}'], text=True).strip() == HEAD
assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
review = json.loads((OLD / 'controller-review.json').read_text())
assert review['accepted']
records = {}
for name in ('run_native_control.py', 'run_shader_clock_control.ps1',
             'emergency_native_cleanup.ps1', 'run_native_retirement.py',
             'verify_native_retirement.ps1', 'prepare_native_stage.py'):
    source = OLD / name
    assert seal(source) == review['controller_seals'][name]
    text = source.read_text()
    if name == 'run_native_control.py':
        assert text.count("HERE/'source/source-review-final.json'") == 1
        text = text.replace("HERE/'source/source-review-final.json'", "HERE/'host-source-review.json'")
    if name == 'prepare_native_stage.py':
        text = text.replace('8f1118dd56d10cdd0f932613427b68409a1f77f6', HEAD)
        text = text.replace('native-8f1118d-shader-clock-stream', 'native-24a8e14-rhs-layer-clock')
        text = text.replace("path = ROOT / 'attestations/native-vulkan/2053ef7/native-shader-clock-admission/files/phase' / name",
                            "path = (WORK / 'phase/retained_transport.spv' if name == 'retained_transport.spv' else ROOT / 'attestations/native-vulkan/2053ef7/native-shader-clock-admission/files/phase/retained_endpoint.spv')")
        text = text.replace("'out/native-shader-clock-execution/source/source-review-final.json', 'out/native-shader-clock-execution/prepare_native_stage.py'",
                            "'out/native-rhs-layer-clock-admission/host-source-review.json', 'out/native-rhs-layer-clock-admission/prepare_native_stage.py', 'out/native-rhs-layer-clock-admission/phase/retained_transport.spv', 'out/native-rhs-layer-clock-admission/compiled-review.json', 'out/native-rhs-layer-clock-admission/compiler.json', 'out/native-rhs-layer-clock-admission/preparation.json'")
        text = text.replace("len(inputs) - 5", "len(inputs) - 9")
        text = text.replace('uncalibrated native-default shader-clock sidecars;', 'uncalibrated native-default shader-clock phase/layer sidecars and exact input snapshots;')
    target = WORK / name
    assert not target.exists()
    target.write_text(text)
    records[name] = {'basis': str(source.relative_to(ROOT)), 'basis_seal': seal(source),
                     'current_seal': seal(target), 'byte_identical': source.read_bytes() == target.read_bytes()}
(WORK / 'native-source-preparation.json').write_text(json.dumps({'source_revision': HEAD, 'records': records}, indent=2) + '\n')
print(json.dumps({'source_revision': HEAD, 'sources_prepared': list(records)}))
