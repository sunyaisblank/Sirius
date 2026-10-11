"""Open the sealed saved trace and select the observed owned EVENTS tab."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
WORK = Path(__file__).resolve().parent
def replace(s, old, new):
    assert s.count(old) == 1, (old, s.count(old))
    return s.replace(old, new)
def save(name, s):
    p = WORK / name
    assert not p.exists()
    p.write_text(s)
ps = (WORK / 'profiler_trace_stable.ps1').read_text().replace(
    "'profiler-trace-stable'", "'profiler-trace-events'")
ps = replace(ps, "    if ($readyRoots -ne 1) {Start-Sleep -Milliseconds 200;continue}", """    if ($readyRoots -ne 1) {Start-Sleep -Milliseconds 200;continue}
    if ($rows.Count -ne 1 -or $process.HasExited) {throw 'Navigation requires exactly one live owned root'}
    $root=[System.Windows.Automation.AutomationElement]::FromHandle([IntPtr]$rows[0].handle)
    if ($root.Current.ProcessId -ne $process.Id) {throw 'Navigation root owner mismatch'}
    $tabName=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::NameProperty,'EVENTS')
    $tabType=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::ControlTypeProperty,[System.Windows.Automation.ControlType]::TabItem)
    $tabCondition=New-Object System.Windows.Automation.AndCondition($tabName,$tabType)
    $tabs=$root.FindAll([System.Windows.Automation.TreeScope]::Descendants,$tabCondition)
    if ($tabs.Count -ne 1) {throw 'Expected unique observed EVENTS tab'}
    $tab=$tabs[0];$c=$tab.Current
    if ($process.HasExited -or $c.ProcessId -ne $process.Id -or !$c.IsEnabled -or $c.IsOffscreen) {throw 'EVENTS tab owner/lifetime/state mismatch'}
    $selection=$null
    if (!$tab.TryGetCurrentPattern([System.Windows.Automation.SelectionItemPattern]::Pattern,[ref]$selection)) {throw 'EVENTS selection pattern unavailable'}
    $record.navigation_requested='EVENTS SelectionItemPattern.Select on exact owned root'
    $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
    if ($process.HasExited -or $root.Current.ProcessId -ne $process.Id -or $tab.Current.ProcessId -ne $process.Id) {throw 'Navigation ownership changed before select'}
    $selection.Select()
    Start-Sleep -Milliseconds 500
    if ($process.HasExited -or $root.Current.ProcessId -ne $process.Id -or $tab.Current.ProcessId -ne $process.Id -or !$selection.Current.IsSelected) {throw 'EVENTS selection did not remain owned/selected'}
    $record.navigation_completed=$true
""")
ps = replace(ps, "if (!$record.loaded_profile_observed -or $record.status", "if (!$record.navigation_completed -or !$record.loaded_profile_observed -or $record.status")
ps = replace(ps, "    if (!$record.loaded_profile_observed) {Start-Sleep -Milliseconds 200;continue}", "    if (!$record.loaded_profile_observed) {throw 'Loaded profile readiness lost after EVENTS selection'}")
ps = replace(ps, "Owned sealed-profile GUI observation only;", "Owned sealed-profile EVENTS tab selection and GUI observation only;")
save('profiler_trace_events.ps1', ps)
b = (WORK / 'run_profiler_trace_stable.py').read_text().replace(
    'profiler-trace-stable', 'profiler-trace-events').replace(
    'profiler_trace_stable.ps1', 'profiler_trace_events.ps1').replace(
    'trace-stable-source-review.json', 'trace-events-source-review.json').replace(
    'prepare_trace_stable.py', 'prepare_trace_events.py')
b = replace(b, "o['terminal_returncode']==2**32-1 and", "o['navigation_completed'] and o['terminal_returncode']==2**32-1 and")
b = replace(b, "Owned GUI window observation only.", "Owned saved-profile EVENTS tab selection and window observation only.")
b = replace(b, "WORK/'profiler-trace-ui-bridge.json'", "WORK/'profiler-trace-ui-bridge.json',WORK/'profiler-trace-stable-bridge.json',WORK/'profiler-trace-stable/automation.json'")
save('run_profiler_trace_events.py', b)
