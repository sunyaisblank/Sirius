"""Use the observed dispatch InvokePattern and inspect its actual effect."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
WORK = Path(__file__).resolve().parent
def replace(s, old, new):
    assert s.count(old) == 1, (old, s.count(old))
    return s.replace(old, new)
def save(n, s):
    p = WORK / n
    assert not p.exists()
    p.write_text(s)
ps = (WORK / 'profiler_trace_selection.ps1').read_text().replace(
    'profiler-trace-selection', 'profiler-trace-invoke')
ps = replace(ps, '$event.TryGetCurrentPattern([System.Windows.Automation.SelectionItemPattern]::Pattern', '$event.TryGetCurrentPattern([System.Windows.Automation.InvokePattern]::Pattern')
ps = replace(ps, '     $selection.Select();Start-Sleep -Milliseconds 500', '     $selection.Invoke();Start-Sleep -Milliseconds 500')
ps = replace(ps, ';original_is_selected=$selection.Current.IsSelected', ";requested_action='InvokePattern.Invoke'")
ps = replace(ps, 'Owned saved-profile tabs and one dispatch-selection request/result observation only; requested selection is not assumed successful;', 'Owned saved-profile tabs and one observed dispatch InvokePattern request/result observation only; its effect is not assumed;')
save('profiler_trace_invoke.ps1', ps)
b = (WORK / 'run_profiler_trace_selection.py').read_text().replace(
    'profiler-trace-selection', 'profiler-trace-invoke').replace(
    'profiler_trace_selection.ps1', 'profiler_trace_invoke.ps1').replace(
    'trace-selection-source-review.json', 'trace-invoke-source-review.json').replace(
    'prepare_trace_selection.py', 'prepare_trace_invoke.py')
b = replace(b, 'Owned saved-profile tabs and one dispatch-selection request/result observation only; no successful selection assumed.', 'Owned saved-profile tabs and one dispatch InvokePattern request/result observation only; no effect assumed.')
b = replace(b, "WORK/'profiler-trace-details/owner.json'", "WORK/'profiler-trace-details/owner.json',WORK/'profiler-trace-selection-bridge.json',WORK/'profiler-trace-selection/owner.json'")
save('run_profiler_trace_invoke.py', b)
