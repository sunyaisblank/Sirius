"""Select both observed dispatch rows and preserve their independent details."""
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
ps = (WORK / 'profiler_trace_timing.ps1').read_text().replace(
    'profiler-trace-timing', 'profiler-trace-details')
start = ps.index("    if ($process.HasExited) {throw 'Owned profiler exited during window observation'}")
end = ps.index('    $record.intentional_owned_stop=$true', start)
inventory = ps[start:end]
inventory = replace(inventory, "(Join-Path $Output 'automation.json')", "(Join-Path $Output ($label+'-automation.json'))")
inventory = inventory.replace('$record.window_observation=$rows', '$record.current_window_observation=$rows')
select = """    $record.event_observations=@()
    foreach ($eventName in @('4 vkCmdDispatch(3, 1, 1)','10 vkCmdDispatch(6, 1, 1)')) {
     $label=if ($eventName -eq '4 vkCmdDispatch(3, 1, 1)') {'transport'} else {'endpoint'}
     if ($process.HasExited -or $root.Current.ProcessId -ne $process.Id) {throw 'Event selection root owner/lifetime mismatch'}
     $eventNameCondition=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::NameProperty,$eventName)
     $eventTypeCondition=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::ControlTypeProperty,[System.Windows.Automation.ControlType]::TreeItem)
     $eventCondition=New-Object System.Windows.Automation.AndCondition($eventNameCondition,$eventTypeCondition)
     $events=$root.FindAll([System.Windows.Automation.TreeScope]::Descendants,$eventCondition)
     if ($events.Count -ne 1) {throw 'Observed dispatch must be unique'}
     $event=$events[0];$c=$event.Current
     if ($process.HasExited -or $c.ProcessId -ne $process.Id -or !$c.IsEnabled -or $c.IsOffscreen) {throw 'Dispatch selection owner/lifetime/state mismatch'}
     $selection=$null
     if (!$event.TryGetCurrentPattern([System.Windows.Automation.SelectionItemPattern]::Pattern,[ref]$selection)) {throw 'Dispatch selection pattern unavailable'}
     $record.event_selection_requested=$eventName
     $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
     if ($process.HasExited -or $event.Current.ProcessId -ne $process.Id -or $root.Current.ProcessId -ne $process.Id) {throw 'Dispatch selection ownership changed'}
     $selection.Select();Start-Sleep -Milliseconds 500
     if ($process.HasExited -or $event.Current.ProcessId -ne $process.Id -or $root.Current.ProcessId -ne $process.Id -or !$selection.Current.IsSelected) {throw 'Dispatch selection did not remain owned/selected'}
"""
finish = """    foreach ($row in $rows) {
     if ($row.print_window_succeeded -and $row.image_path) {
      $renamed=Join-Path $Output ($label+'-'+[IO.Path]::GetFileName($row.image_path))
      if (Test-Path -LiteralPath $renamed) {throw 'One-shot event image already exists'}
      Move-Item -LiteralPath $row.image_path -Destination $renamed
      $row.image_path=$renamed
     }
    }
    $record.event_observations+=@([pscustomobject]@{event_name=$eventName;label=$label;owned_ui_elements=$automation.Count;automation_path=(Join-Path $Output ($label+'-automation.json'));window_observation=$rows})
    }
    $record.window_observation=@($record.event_observations | ForEach-Object {$_.window_observation})
"""
ps = ps[:start] + select + inventory + finish + ps[end:]
ps = replace(ps, 'Owned sealed-profile EVENTS and Event timing tab selection and GUI observation only;', 'Owned sealed-profile EVENTS/Event timing and two observed dispatch selections with GUI observation only;')
ps = replace(ps, "if (!$record.navigation_completed -or", "if ($record.event_observations.Count -ne 2 -or !$record.navigation_completed -or")
save('profiler_trace_details.ps1', ps)
b = (WORK / 'run_profiler_trace_timing.py').read_text().replace(
    'profiler-trace-timing', 'profiler-trace-details').replace(
    'profiler_trace_timing.ps1', 'profiler_trace_details.ps1').replace(
    'trace-timing-source-review.json', 'trace-details-source-review.json').replace(
    'prepare_trace_timing.py', 'prepare_trace_details.py')
b = replace(b, "o['navigation_history']==", "[v['event_name'] for v in o['event_observations']]==['4 vkCmdDispatch(3, 1, 1)','10 vkCmdDispatch(6, 1, 1)'] and o['navigation_history']==")
b = replace(b, 'Owned saved-profile EVENTS and Event timing selection and window observation only.', 'Owned saved-profile EVENTS/Event timing and two observed dispatch selections with window observation only.')
b = replace(b, "WORK/'profiler-trace-events/automation.json'", "WORK/'profiler-trace-events/automation.json',WORK/'profiler-trace-timing-bridge.json',WORK/'profiler-trace-timing/automation.json'")
save('run_profiler_trace_details.py', b)
