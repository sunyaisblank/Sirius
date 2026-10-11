"""Move only acquisition to the first legacy T/E pair; keep the science unchanged."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import json
WORK=Path(__file__).resolve().parent
def replace(s,old,new):
 assert s.count(old)==1,(old,s.count(old));return s.replace(old,new)
def save(n,s):
 p=WORK/n;assert not p.exists();p.write_text(s)
plan=json.loads((WORK/'capture-plan.json').read_text())
plan.update(dispatch_spec='dispatch:9:2',expected_dispatch_dimensions=[[3,1,1],[6,1,1]],expected_marker_indices=[0,1],expected_phase='legacy_repeat0',automatic_profiler_peak_clocks=True,clock_scope='RDP tracing automatically requests peak clocks and restoration. No explicit clock or driver experiment configuration; actual frequencies/restoration not independently measured. No unprofiled benefit comparison.',scope='One changed early legacy-range acquisition on the unchanged original numerical control. More natural remaining device lifetime is a hypothesis, not a proven failure cause or guaranteed dump margin. Actual trace must confirm events/pipeline identity and sample availability.')
save('capture-early-plan.json',json.dumps(plan,indent=2)+'\n')
ps=(WORK/'run_console_capture_control.ps1').read_text()
ps=replace(ps,'--rgp-auto-capture=dispatch:43:2','--rgp-auto-capture=dispatch:9:2')
ps=replace(ps,'original-modern-te.rgp','original-legacy-te.rgp')
ps=replace(ps,'with dispatch:43:2 instruction/counter trace.','with dispatch:9:2 legacy instruction/counter trace.')
ps=replace(ps,'No clocks, driver experiments, shader instrumentation, performance benefit or full qualification claim.','No explicit clock or driver experiment configuration. RDP automatically requests peak clocks and restoration; actual frequencies/restored state unmeasured. No shader instrumentation, performance benefit or full qualification claim.')
save('run_early_capture_control.ps1',ps)
reader=(WORK/'verify_console_capture_control.py').read_text()
reader=replace(reader,"WORK/'capture-plan.json'","WORK/'capture-early-plan.json'")
assert reader.count('original-modern-te.rgp')==2;reader=reader.replace('original-modern-te.rgp','original-legacy-te.rgp')
reader=replace(reader,'/native-capture-console/','/native-capture-early-console/')
save('verify_early_capture_control.py',reader)
bridge=(WORK/'run_native_console_capture.py').read_text()
bridge=bridge.replace('native-capture-console','native-capture-early-console')
bridge=bridge.replace('native-capture-early-console-source-review.json','native-capture-early-source-review.json')
bridge=bridge.replace('run_console_capture_control.ps1','run_early_capture_control.ps1')
bridge=bridge.replace('verify_console_capture_control','verify_early_capture_control')
bridge=bridge.replace("WORK/'capture-plan.json'","WORK/'capture-early-plan.json'")
bridge=replace(bridge,'original-modern-te.rgp','original-legacy-te.rgp')
bridge=replace(bridge,"paths += [WORK/'owned_console_capture.cs'","paths += [WORK/'prepare_early_capture.py',WORK/'verify_failed_capture_result.py',WORK/'native-capture-console-result.json',WORK/'owned_console_capture.cs'")
bridge=replace(bridge,'complete readback and stream verification;', 'early legacy range; complete readback and stream verification; automatic RDP peak/restore clock requests are reported tool effects, not ordinary performance evidence;')
save('run_native_early_capture.py',bridge)
