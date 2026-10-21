from pathlib import Path
import json,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;label,folder=sys.argv[1:];a=json.loads((W/'actions.json').read_text())[label]
for name in ['candidate-default','candidate-poison']:
    v=json.loads((R/'out/native-endpoint-phase-sharing/reused/prior-proof/lavapipe-ieee32-validation'/name/'numerical-and-boundary-result.json').read_text());assert v['numerical_pass'] and v['all_observed_callback_states_restored']
assert json.loads((R/'out/native-endpoint-phase-sharing/ordinary/ordinary-code-review.json').read_text())['structural_premise_observed']
for previous in a['prerequisites']:
    v=json.loads((W/previous/'owner.json').read_text());assert v['passed'] and v['remaining_owned_processes']==[]
if label.startswith('control-') and label!='control-8':
    v=json.loads((R/'out/native-endpoint-phase-sharing/readers/benefit.json').read_text());assert v['audit_passed'] and v['all22_numeric_pass']
subprocess.run([str(R/'bin/linux-gcc/native-endpoint-phase-sharing/matched/audited_backend_tests'),'--gtest_filter='+a['case'],'--gtest_color=no'],check=True)
