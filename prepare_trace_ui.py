"""Observe the sealed early trace in an owned RGP window, without recapturing."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
WORK=Path(__file__).resolve().parent
def replace(s,old,new):
 assert s.count(old)==1,(old,s.count(old));return s.replace(old,new)
def save(n,s):
 p=WORK/n;assert not p.exists();p.write_text(s)
ps=(WORK/'profiler_ui.ps1').read_text()
ps=replace(ps,"[ValidateSet('profiler-ui')]","[ValidateSet('profiler-trace-ui')]")
ps=replace(ps,'[Parameter(Mandatory=$true)][string]$LiveRevision)','[Parameter(Mandatory=$true)][string]$LiveRevision,[Parameter(Mandatory=$true)][string]$Trace,[Parameter(Mandatory=$true)][string]$TraceHash)')
ps=replace(ps,'Add-Type -ReferencedAssemblies System.Drawing,System,System.Core', "Add-Type -AssemblyName UIAutomationClient,UIAutomationTypes\nif (!(Test-Path -LiteralPath $Trace) -or (Get-Item -LiteralPath $Trace).Length -le 0 -or (Get-Item -LiteralPath $Trace).Length -gt 268435456 -or (Get-FileHash -LiteralPath $Trace -Algorithm SHA256).Hash.ToLowerInvariant() -ne $TraceHash) {throw 'Pinned trace input mismatch'}\nAdd-Type -ReferencedAssemblies System.Drawing,System,System.Core")
ps=replace(ps,"arguments=@();intentional_owned_stop", "arguments=@($Trace);trace_sha256=$TraceHash;trace_path=$Trace;intentional_owned_stop")
ps=replace(ps,'Owned profiler GUI feasibility observation only.', 'Owned sealed-profile GUI observation only; positional file opening is a bounded probe of the documented executable-drop route.')
ps=replace(ps,'-FilePath $exe -WorkingDirectory $Tools', '-FilePath $exe -ArgumentList @((\'"\'+$Trace+\'"\')) -WorkingDirectory $Tools')
ps=replace(ps,'$watch.Elapsed.TotalSeconds -ge 3','$watch.Elapsed.TotalSeconds -ge 5')
ps=replace(ps,'    $record.window_observation=$rows', """    $record.window_observation=$rows
    $automation=New-Object 'System.Collections.Generic.List[object]'
    foreach ($row in $rows) {
     if ($process.HasExited) {throw 'Owned profiler exited before UI observation'}
     $root=[System.Windows.Automation.AutomationElement]::FromHandle([IntPtr]$row.handle)
     if ($root.Current.ProcessId -ne $process.Id) {throw 'UI root owner mismatch'}
     $condition=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::ProcessIdProperty,[int]$process.Id)
     $elements=$root.FindAll([System.Windows.Automation.TreeScope]::Descendants,$condition)
     if ($automation.Count+$elements.Count -gt 5000) {throw 'Owned UI inventory exceeds diagnostic bound'}
     foreach ($element in $elements) {
      $c=$element.Current
      if ($c.ProcessId -ne $process.Id -or $process.HasExited) {throw 'UI child owner/lifetime mismatch'}
      $value=$null;$pattern=$null
      if ($element.TryGetCurrentPattern([System.Windows.Automation.ValuePattern]::Pattern,[ref]$pattern)) {$value=$pattern.Current.Value}
      if ($c.Name.Length -gt 4096 -or ($null -ne $value -and $value.Length -gt 4096)) {throw 'Owned UI text exceeds diagnostic bound'}
      $automation.Add([pscustomobject]@{name=$c.Name;value=$value;control=$c.ControlType.ProgrammaticName;automation_id=$c.AutomationId;class=$c.ClassName;enabled=$c.IsEnabled;offscreen=$c.IsOffscreen;rect=@($c.BoundingRectangle.X,$c.BoundingRectangle.Y,$c.BoundingRectangle.Width,$c.BoundingRectangle.Height)})
     }
    }
    $automation.ToArray() | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'automation.json') -Encoding UTF8
    $record.owned_ui_elements=$automation.Count""")
ps=replace(ps,' $record.tool_unchanged=$false', " $record.trace_unchanged=$false\n try {$record.trace_unchanged=((Get-FileHash -LiteralPath $Trace -Algorithm SHA256).Hash.ToLowerInvariant() -eq $TraceHash)} catch {$errors.Add('trace-posthash: '+$_.Exception.Message)}\n $record.tool_unchanged=$false")
ps=replace(ps,' if ($errors.Count -or !$record.owned_process_absent', ' if (!$record.trace_unchanged -or $errors.Count -or !$record.owned_process_absent')
save('profiler_trace_ui.ps1',ps)
b=(WORK/'run_profiler_ui.py').read_text()
b=b.replace("'profiler-ui'","'profiler-trace-ui'")
b=b.replace('gui-inspection-source-review.json','trace-ui-source-review.json')
b=b.replace('profiler_ui.ps1','profiler_trace_ui.ps1')
b=replace(b,"paths += [ROOT/v['path']", "trace=WORK/'native-capture-early-console/original-legacy-te.rgp';acquisition=read(WORK/'native-capture-early-console-bridge.json');assert acquisition['accepted']\nassert trace.stat().st_size==acquisition['original_control_and_acquired_trace']['trace_bytes'] and sha(trace)==acquisition['original_control_and_acquired_trace']['trace_sha256']\npaths += [trace,WORK/'native-capture-early-console-bridge.json',WORK/'prepare_trace_ui.py']\npaths += [ROOT/v['path']")
b=replace(b,"'-LiveRevision',head]", "'-LiveRevision',head,'-Trace',win(trace),'-TraceHash',sha(trace)]")
b=replace(b,"assert o['test_executable_sha256']==sha(exe)","assert o['test_executable_sha256']==sha(exe) and o['trace_unchanged'] and o['trace_sha256']==sha(trace)")
save('run_profiler_trace_ui.py',b)
