"""Observe both dispatch details and instruction panes of the same sealed trace."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import json
WORK = Path(__file__).resolve().parent
def replace(s, old, new):
    assert s.count(old) == 1, (old, s.count(old))
    return s.replace(old, new)
def save(n, s):
    p = WORK / n
    assert not p.exists()
    p.write_text(s)
ps = (WORK / 'profiler_trace_owned_click.ps1').read_text().replace(
    'profiler-trace-owned-click', 'profiler-trace-samples')
start = ps.index('    $tabName=')
end = ps.index('    }\n    $record.navigation_completed=$true', start)
body = ps[start:end]
function = 'function SelectOwnedTab([string]$tabLabel) {\n' + body + '}\n'
ps = replace(ps, body, '    SelectOwnedTab $tabLabel\n')
function = replace(function, '    $tabName=', '    GuardOwnedView\n    $tabName=')
guard = """function GuardOwnedView {
 $process.Refresh()
 if ($process.HasExited) {throw 'Owned profiler exited before view action'}
 $script:peak=[Math]::Max($script:peak,$process.WorkingSet64)
 if ($watch.Elapsed.TotalSeconds -gt 60 -or $process.WorkingSet64 -gt 1073741824) {$record.outer_stop=$true;throw 'Owned view guard exceeded'}
}
"""
ps = replace(ps, '$process=$null;$watch=', guard+function+'$process=$null;$watch=')
ps = replace(ps, "foreach ($eventName in @('4 vkCmdDispatch(3, 1, 1)'))", "foreach ($eventName in @('4 vkCmdDispatch(3, 1, 1)','10 vkCmdDispatch(6, 1, 1)'))")
ps = replace(ps, "     $label=if ($eventName", "     GuardOwnedView\n     if ($eventName -eq '10 vkCmdDispatch(6, 1, 1)') {SelectOwnedTab 'Event timing'}\n     $baseLabel=if ($eventName")
a = json.loads((WORK/'profiler-trace-owned-click/transport-automation.json').read_text(encoding='utf-8-sig'))
event_ids = [v['automation_id'] for v in a if v['automation_id'].endswith('SingleEventSidePanel.event_name_')]
assert len(event_ids) == 1
before_inventory = ps.index("    if ($process.HasExited) {throw 'Owned profiler exited during window observation'}")
setup = """     if (!$record.selection_result.fresh_is_selected) {throw 'Owned native click did not select dispatch'}
     $identityCondition=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::AutomationIdProperty,'EVENT_ID')
     $identity=$root.FindFirst([System.Windows.Automation.TreeScope]::Descendants,$identityCondition)
     if (!$identity -or $process.HasExited -or $identity.Current.ProcessId -ne $process.Id -or $identity.Current.Name -ne $eventName) {throw 'Selected event detail identity mismatch'}
     foreach ($view in @('event','instructions')) {
      GuardOwnedView
      $label=$baseLabel+'-'+$view
      if ($view -eq 'instructions') {SelectOwnedTab 'Instruction timing'}
      if ($process.HasExited) {throw 'Owned profiler exited before view screenshot'}
      $rows=@([SiriusWindow]::Capture([uint32]$process.Id,$Output))
      $record.last_observation_label=$label;$record.current_window_observation=$rows
      $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
""".replace('EVENT_ID', event_ids[0])
ps = ps[:before_inventory]+setup+ps[before_inventory:]
ps = replace(ps, 'event_name=$eventName;label=$label;owned_ui_elements=', 'event_name=$eventName;view=$view;label=$label;selected_event_confirmed=$true;selection_state=$record.selection_result;owned_ui_elements=')
ps = replace(ps, '    }\n    $record.window_observation=@($record.event_observations', '    }\n    }\n    $record.window_observation=@($record.event_observations')
ps = replace(ps, '$record.event_observations.Count -ne 1', '$record.event_observations.Count -ne 4')
ps = replace(ps, 'Owned saved-profile tabs and one bounded HWND client mouse-message request/result observation only; no global input API or assumed effect;', 'Owned saved-profile tabs and two bounded HWND dispatch selections with independent event and instruction-pane observations; selected event identity required; no global input API or science;')
save('profiler_trace_samples.ps1', ps)
b = (WORK / 'run_profiler_trace_owned_click.py').read_text().replace(
    'profiler-trace-owned-click', 'profiler-trace-samples').replace(
    'profiler_trace_owned_click.ps1', 'profiler_trace_samples.ps1').replace(
    'trace-owned-click-source-review.json', 'trace-samples-source-review.json').replace(
    'prepare_trace_owned_click.py', 'prepare_trace_samples.py')
b = replace(b, "[v['event_name'] for v in o['event_observations']]==['4 vkCmdDispatch(3, 1, 1)']", "[(v['event_name'],v['view']) for v in o['event_observations']]==[('4 vkCmdDispatch(3, 1, 1)','event'),('4 vkCmdDispatch(3, 1, 1)','instructions'),('10 vkCmdDispatch(6, 1, 1)','event'),('10 vkCmdDispatch(6, 1, 1)','instructions')] and all(v['selected_event_confirmed'] and v['selection_state']['fresh_is_selected'] for v in o['event_observations'])")
b = replace(b, "o['navigation_history']==['EVENTS','Event timing']", "o['navigation_history']==['EVENTS','Event timing','Instruction timing','Event timing','Instruction timing']")
b = replace(b, 'Owned saved-profile tabs and one bounded owned-HWND client mouse-message request/result observation only; no global input or effect assumed.', 'Owned saved-profile tabs and two bounded owned-HWND dispatch selections; independent event and instruction-pane observations; no sample interpretation or science assumed.')
b = replace(b, "WORK/'profiler-trace-invoke/owner.json'", "WORK/'profiler-trace-invoke/owner.json',WORK/'profiler-trace-owned-click-bridge.json',WORK/'profiler-trace-owned-click/transport-automation.json'")
save('run_profiler_trace_samples.py', b)
