from pathlib import Path
import json,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;label,folder=sys.argv[1:];folder=Path(folder)
a=json.loads((W/'actions.json').read_text())[label]
assert json.loads((R/'out/portable-checked-fp64-matched/matched-result.json').read_text())['passed']
assert json.loads((R/'out/portable-checked-fp64-jit/ordinary-code-review.json').read_text())['structural_premise_observed']
assert json.loads((R/'out/portable-checked-fp64/run/complete-output-agreement.json').read_text())['passed']
for previous in a['prerequisites']:
 v=json.loads((W/previous/'owner.json').read_text());assert v['passed'] and v['remaining_owned_processes']==[]
subprocess.run([str(R/f'bin/linux-gcc/portable-checked-fp64-matched/{a["arm"]}-tests'),'--gtest_filter=RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults','--gtest_color=no'],check=True)
