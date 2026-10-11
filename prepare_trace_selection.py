"""Observe the actual result of one dispatch-selection request; no success assumption."""
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
ps = (WORK / 'profiler_trace_details.ps1').read_text().replace(
    'profiler-trace-details', 'profiler-trace-selection')
ps = replace(ps, "foreach ($eventName in @('4 vkCmdDispatch(3, 1, 1)','10 vkCmdDispatch(6, 1, 1)'))", "foreach ($eventName in @('4 vkCmdDispatch(3, 1, 1)'))")
ps = replace(ps, "     if ($process.HasExited -or $event.Current.ProcessId -ne $process.Id -or $root.Current.ProcessId -ne $process.Id -or !$selection.Current.IsSelected) {throw 'Dispatch selection did not remain owned/selected'}", """     if ($process.HasExited -or $root.Current.ProcessId -ne $process.Id) {throw 'Owned root/lifetime changed after selection request'}
     $record.selection_result=[ordered]@{original_event_pid=$event.Current.ProcessId;original_is_selected=$selection.Current.IsSelected}
     $fresh=$root.FindAll([System.Windows.Automation.TreeScope]::Descendants,$eventCondition)
     $record.selection_result.fresh_count=$fresh.Count
     if ($fresh.Count -ne 1 -or $process.HasExited -or $fresh[0].Current.ProcessId -ne $process.Id) {throw 'Fresh dispatch reference is not uniquely owned'}
     $freshSelection=$null
     if (!$fresh[0].TryGetCurrentPattern([System.Windows.Automation.SelectionItemPattern]::Pattern,[ref]$freshSelection)) {throw 'Fresh dispatch selection pattern unavailable'}
     $record.selection_result.fresh_is_selected=$freshSelection.Current.IsSelected
     $record.selection_result.available_patterns=@($fresh[0].GetSupportedPatterns() | ForEach-Object {$_.ProgrammaticName})
     if ($process.HasExited -or $root.Current.ProcessId -ne $process.Id -or $fresh[0].Current.ProcessId -ne $process.Id) {throw 'Fresh selection observation ownership changed'}
""")
ps = replace(ps, '$record.event_observations.Count -ne 2', '$record.event_observations.Count -ne 1')
ps = replace(ps, 'Owned sealed-profile EVENTS/Event timing and two observed dispatch selections with GUI observation only;', 'Owned saved-profile tabs and one dispatch-selection request/result observation only; requested selection is not assumed successful;')
save('profiler_trace_selection.ps1', ps)
b = (WORK / 'run_profiler_trace_details.py').read_text().replace(
    'profiler-trace-details', 'profiler-trace-selection').replace(
    'profiler_trace_details.ps1', 'profiler_trace_selection.ps1').replace(
    'trace-details-source-review.json', 'trace-selection-source-review.json').replace(
    'prepare_trace_details.py', 'prepare_trace_selection.py')
b = replace(b, "==['4 vkCmdDispatch(3, 1, 1)','10 vkCmdDispatch(6, 1, 1)']", "==['4 vkCmdDispatch(3, 1, 1)']")
b = replace(b, 'Owned saved-profile EVENTS/Event timing and two observed dispatch selections with window observation only.', 'Owned saved-profile tabs and one dispatch-selection request/result observation only; no successful selection assumed.')
b = replace(b, "WORK/'profiler-trace-timing/automation.json'", "WORK/'profiler-trace-timing/automation.json',WORK/'profiler-trace-details-bridge.json',WORK/'profiler-trace-details/owner.json'")
save('run_profiler_trace_selection.py', b)
