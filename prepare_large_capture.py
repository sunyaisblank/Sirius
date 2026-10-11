"""Acquire the original disabled-modern pair with more waves; change no science."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import json
WORK = Path(__file__).resolve().parent
outputs = {}
def replace(s, old, new):
    assert s.count(old) == 1, (old, s.count(old))
    return s.replace(old, new)
def save(name, text):
    path = WORK / name
    assert not path.exists() and name not in outputs
    outputs[name] = text
plan = json.loads((WORK / 'capture-early-plan.json').read_text())
plan.update(
    dispatch_spec='dispatch:32:2',
    expected_dispatch_dimensions=[[12, 1, 1], [24, 1, 1]],
    expected_marker_indices=[],
    expected_phase='disabled_modern',
    scope='One changed acquisition range on the unchanged original numerical control. The disabled-modern Transport/Endpoint pair has more waves than legacy dispatches 9/10, whose actual instruction panes yielded no usable samples. Range 32/33 has no measured timestamp-marker pair; all 51 original measured markers must still pass. More waves and remaining natural device lifetime may help acquisition, but neither sample availability nor dump margin is guaranteed. Actual trace event/pipeline identity and samples require separate verification.')
save('capture-large-plan.json', json.dumps(plan, indent=2) + '\n')
ps = (WORK / 'run_early_capture_control.ps1').read_text()
ps = replace(ps, '--rgp-auto-capture=dispatch:9:2', '--rgp-auto-capture=dispatch:32:2')
ps = replace(ps, 'original-legacy-te.rgp', 'original-disabled-modern-te.rgp')
ps = replace(ps, 'with dispatch:9:2 legacy instruction/counter trace.',
             'with dispatch:32:2 disabled-modern instruction/counter trace. No measured marker pair belongs to this disabled range.')
save('run_large_capture_control.ps1', ps)
reader = (WORK / 'verify_early_capture_control.py').read_text()
reader = replace(reader, "WORK/'capture-early-plan.json'", "WORK/'capture-large-plan.json'")
assert reader.count('original-legacy-te.rgp') == 2
reader = reader.replace('original-legacy-te.rgp', 'original-disabled-modern-te.rgp')
reader = replace(reader, '/native-capture-early-console/', '/native-capture-large-console/')
save('verify_large_capture_control.py', reader)
bridge = (WORK / 'run_native_early_capture.py').read_text()
bridge = bridge.replace('native-capture-early-console', 'native-capture-large-console')
bridge = bridge.replace('native-capture-early-source-review.json', 'native-capture-large-source-review.json')
bridge = bridge.replace('run_early_capture_control.ps1', 'run_large_capture_control.ps1')
bridge = bridge.replace('verify_early_capture_control', 'verify_large_capture_control')
assert bridge.count("WORK/'capture-early-plan.json'") == 2
bridge = bridge.replace("WORK/'capture-early-plan.json'", "WORK/'capture-large-plan.json'")
bridge = replace(bridge, 'original-legacy-te.rgp', 'original-disabled-modern-te.rgp')
bridge = replace(bridge, "WORK/'prepare_early_capture.py'",
                 "WORK/'prepare_large_capture.py',WORK/'native-capture-early-console-bridge.json',WORK/'profiler-trace-samples-bridge.json',WORK/'profiler-trace-samples/transport-instructions-automation.json',WORK/'profiler-trace-samples/endpoint-instructions-automation.json'")
bridge = replace(bridge, 'early legacy range;',
                 'larger disabled-modern range with no associated measured marker pair;')
save('run_native_large_capture.py', bridge)
for name, text in outputs.items():
    (WORK / name).write_text(text)
