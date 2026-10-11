"""Capture all original legacy-repeat0 dispatches without changing the test."""
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
    assert not (WORK / name).exists() and name not in outputs
    outputs[name] = text
plan = json.loads((WORK / 'capture-early-plan.json').read_text())
plan.update(dispatch_spec='dispatch:9:7',
            expected_dispatch_dimensions=[[3,1,1],[6,1,1],[3,1,1],[6,1,1],
                                          [3,1,1],[6,1,1],[3,1,1]],
            expected_marker_indices=list(range(6)),
            expected_phase='legacy_repeat0',
            scope='One changed coverage range on the unchanged original numerical control. Capture the seven physical dispatches of legacy repeat0: three Transport/Endpoint pairs and final Dense. Six measured marker bundles cover this repeat; marker5 combines Endpoint and Dense. Per-event cohorts remain small; repeated opportunities may still yield no instruction samples. More trace data may reduce dump margin. The last acquired trace has only90.02268ms from marker5 end through marker50 end; no substantial remaining lifetime is established. More natural remaining device lifetime is a hypothesis, not an established failure cause. Range/dimensions alone cannot prove trace-to-marker identity. Actual event/pipeline/code joins and usable samples require separate verification; no ordinary performance or full qualification claim.')
save('capture-interval-plan.json', json.dumps(plan, indent=2) + '\n')
ps = (WORK / 'run_early_capture_control.ps1').read_text()
ps = replace(ps, '--rgp-auto-capture=dispatch:9:2', '--rgp-auto-capture=dispatch:9:7')
ps = replace(ps, "'--rgp-render-op-count=2'", "'--rgp-render-op-count=7'")
ps = replace(ps, 'original-legacy-te.rgp', 'original-legacy-interval.rgp')
ps = replace(ps, 'with dispatch:9:2 legacy instruction/counter trace.',
             'with dispatch:9:7 complete legacy-repeat0 instruction/counter trace. Small per-event cohorts may still have no analyzed samples.')
save('run_interval_capture_control.ps1', ps)
reader = (WORK / 'verify_early_capture_control.py').read_text()
reader = replace(reader, "WORK/'capture-early-plan.json'", "WORK/'capture-interval-plan.json'")
assert reader.count('original-legacy-te.rgp') == 2
reader = reader.replace('original-legacy-te.rgp', 'original-legacy-interval.rgp')
reader = replace(reader, '/native-capture-early-console/', '/native-capture-interval-console/')
save('verify_interval_capture_control.py', reader)
bridge = (WORK / 'run_native_early_capture.py').read_text()
bridge = bridge.replace('native-capture-early-console', 'native-capture-interval-console')
bridge = bridge.replace('native-capture-early-source-review.json', 'native-capture-interval-source-review.json')
bridge = bridge.replace('run_early_capture_control.ps1', 'run_interval_capture_control.ps1')
bridge = bridge.replace('verify_early_capture_control', 'verify_interval_capture_control')
assert bridge.count("WORK/'capture-early-plan.json'") == 2
bridge = bridge.replace("WORK/'capture-early-plan.json'", "WORK/'capture-interval-plan.json'")
bridge = replace(bridge, 'original-legacy-te.rgp', 'original-legacy-interval.rgp')
bridge = replace(bridge, "WORK/'prepare_early_capture.py'",
                 "WORK/'prepare_interval_capture.py',WORK/'trace-samples-result.json',WORK/'native-capture-large-console-result.json'")
bridge = replace(bridge, 'early legacy range;',
                 'seven-dispatch legacy-repeat0 coverage range; per-event sample availability unproved;')
save('run_native_interval_capture.py', bridge)
for name, text in outputs.items():
    (WORK / name).write_text(text)
