"""Observe one pair at a time from the sealed seven-dispatch original trace."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
WORK = Path(__file__).resolve().parent
def replace(s, old, new):
    assert s.count(old) == 1, (old, s.count(old))
    return s.replace(old, new)
ps = (WORK / 'profiler_trace_samples.ps1').read_text()
ps = replace(ps, "[ValidateSet('profiler-trace-samples')]",
             "[ValidateSet('profiler-interval-pair1','profiler-interval-pair2','profiler-interval-pair3')]")
ps = replace(ps, '[string]$TraceHash)',
             '[string]$TraceHash,[Parameter(Mandatory=$true)][ValidateRange(1,3)][int]$Pair)')
ps = replace(ps, "$ErrorActionPreference='Stop'", "$ErrorActionPreference='Stop'\nif ($Action -ne ('profiler-interval-pair'+$Pair)) {throw 'Bounded pair/action mismatch'}")
discovery = """    GuardOwnedView
    $treeCondition=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::ControlTypeProperty,[System.Windows.Automation.ControlType]::TreeItem)
    $treeItems=$root.FindAll([System.Windows.Automation.TreeScope]::Descendants,$treeCondition)
    if ($treeItems.Count -gt 256) {throw 'Dispatch discovery exceeds bound'}
    $dispatches=@()
    foreach ($item in $treeItems) {
     GuardOwnedView
     $current=$item.Current
     if ($current.ProcessId -ne $process.Id -or $root.Current.ProcessId -ne $process.Id) {throw 'Dispatch discovery owner mismatch'}
     if ($current.Name -match '^([0-9]+) vkCmdDispatch\\(([0-9]+), 1, 1\\)$') {
      if (!$current.IsEnabled -or $current.IsOffscreen) {throw 'Dispatch discovery is not available'}
      $dispatches+=@([pscustomobject]@{name=$current.Name;event_index=[int]$Matches[1];groups_x=[int]$Matches[2]})
     }
    }
    $dispatches=@($dispatches | Sort-Object event_index)
    if ($dispatches.Count -ne 7 -or @($dispatches.event_index | Sort-Object -Unique).Count -ne 7 -or ($dispatches.groups_x -join ',') -ne '3,6,3,6,3,6,3') {throw 'Actual seven-dispatch trace dimensions/order mismatch'}
    $record.observed_dispatches=$dispatches;$record.pair=$Pair
    $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
    $chosenEvents=@($dispatches[2*($Pair-1)],$dispatches[2*($Pair-1)+1])
    $record.event_observations=@()
    foreach ($chosen in $chosenEvents) {
     $eventName=$chosen.name
"""
ps = replace(ps, "    $record.event_observations=@()\n    foreach ($eventName in @('4 vkCmdDispatch(3, 1, 1)','10 vkCmdDispatch(6, 1, 1)')) {\n", discovery)
ps = replace(ps, "if ($eventName -eq '10 vkCmdDispatch(6, 1, 1)')", "if ($eventName -eq $chosenEvents[1].name)")
ps = replace(ps, "if ($eventName -eq '4 vkCmdDispatch(3, 1, 1)')", "if ($eventName -eq $chosenEvents[0].name)")
ps = replace(ps, 'selected event identity required;', 'actual seven-dispatch dimension/order discovery and selected event identity required;')
bridge = (WORK / 'run_profiler_trace_samples.py').read_text()
bridge = replace(bridge, "('profiler-trace-samples',)", "('profiler-interval-pair1','profiler-interval-pair2','profiler-interval-pair3')")
bridge = replace(bridge, 'action=sys.argv[1];head=', 'action=sys.argv[1];pair=int(action[-1]);head=')
bridge = bridge.replace('profiler_trace_samples.ps1', 'profiler_interval_pair.ps1')
bridge = bridge.replace('trace-samples-source-review.json', 'interval-observer-source-review.json')
bridge = bridge.replace('prepare_trace_samples.py', 'prepare_interval_observer.py')
bridge = bridge.replace('native-capture-early-console', 'native-capture-interval-console')
bridge = replace(bridge, 'original-legacy-te.rgp', 'original-legacy-interval.rgp')
bridge = replace(bridge, "'-TraceHash',sha(trace)]", "'-TraceHash',sha(trace),'-Pair',str(pair)]")
old = "[(v['event_name'],v['view']) for v in o['event_observations']]==[('4 vkCmdDispatch(3, 1, 1)','event'),('4 vkCmdDispatch(3, 1, 1)','instructions'),('10 vkCmdDispatch(6, 1, 1)','event'),('10 vkCmdDispatch(6, 1, 1)','instructions')]"
new = "o['pair']==pair and len(o['observed_dispatches'])==7 and [v['groups_x'] for v in o['observed_dispatches']]==[3,6,3,6,3,6,3] and len({v['event_index'] for v in o['observed_dispatches']})==7 and [(v['event_name'],v['view']) for v in o['event_observations']]==[(o['observed_dispatches'][i]['name'],view) for i in (2*(pair-1),2*(pair-1)+1) for view in ('event','instructions')]"
bridge = replace(bridge, old, new)
bridge = replace(bridge, 'independent event and instruction-pane observations;', 'actual seven-dispatch dimensions/order recorded; one of three bounded pairs with independent event and instruction-pane observations;')
for name, text in [('profiler_interval_pair.ps1', ps), ('run_profiler_interval_pair.py', bridge)]:
    assert not (WORK / name).exists()
for name, text in [('profiler_interval_pair.ps1', ps), ('run_profiler_interval_pair.py', bridge)]:
    (WORK / name).write_text(text)
