"""Discover current saved-profile events and the observed counter-view controls."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import ast,hashlib,json,subprocess
WORK=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def save(n,s):
 p=WORK/n;assert not p.exists()
 if n.endswith('.py'):ast.parse(s,filename=n)
 p.write_text(s)
def replace(s,o,n):
 assert s.count(o)==1,(o[:80],s.count(o));return s.replace(o,n)
prior=json.loads((WORK/'profiler-counter-inventory-bridge.json').read_text());assert prior['accepted']
u=json.loads((WORK/'profiler-counter-inventory/automation.json').read_text(encoding='utf-8-sig'))
button=[v for v in u if v['name']=='Views' and v['control']=='ControlType.Button'];assert len(button)==1
button=button[0];assert button['enabled'] and not button['offscreen'] and button['class']=='ArrowIconComboBox'
old=subprocess.check_output(['git','show','8b2e8db7efec4bed5ec1de0192b129139545f77a:profiler_interval_second_pair.ps1']).decode()
ps=(WORK/'profiler_counter_inventory.ps1').read_text().replace('profiler-counter-inventory','profiler-counter-controls')
start=old.index(' [StructLayout(LayoutKind.Sequential)] private struct POINT')
end=old.index(' public static SiriusWindowInfo[] Capture',start)
ps=replace(ps,' public static SiriusWindowInfo[] Capture',old[start:end]+' public static SiriusWindowInfo[] Capture')
start=old.index('function GuardOwnedView');end=old.index('$process=$null;',start)
functions=old[start:end]
observe=r'''
function ObserveOwnedView([string]$label) {
 GuardOwnedView
 $rows=@([SiriusWindow]::Capture([uint32]$process.Id,$Output))
 if (!$rows.Count -or $rows.Count -gt 8) {throw 'Owned view window count outside bound'}
 $automation=New-Object 'System.Collections.Generic.List[object]'
 foreach ($row in $rows) {
  GuardOwnedView
  $viewRoot=[System.Windows.Automation.AutomationElement]::FromHandle([IntPtr]$row.handle)
  if ($viewRoot.Current.ProcessId -ne $process.Id) {throw 'View root owner mismatch'}
  $condition=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::ProcessIdProperty,[int]$process.Id)
  $elements=$viewRoot.FindAll([System.Windows.Automation.TreeScope]::Descendants,$condition)
  if ($automation.Count+$elements.Count -gt 5000) {throw 'Owned UI inventory exceeds diagnostic bound'}
  foreach ($element in $elements) {
   GuardOwnedView
   $c=$element.Current
   if ($c.ProcessId -eq 0) {$record.ui_unavailable_pid_zero++;continue}
   if ($c.ProcessId -ne $process.Id) {throw 'UI child foreign owner mismatch'}
   $value=$null;$pattern=$null;$range=$null;$scroll=$null
   if ($element.TryGetCurrentPattern([System.Windows.Automation.ValuePattern]::Pattern,[ref]$pattern)) {$value=$pattern.Current.Value}
   $pattern=$null
   if ($element.TryGetCurrentPattern([System.Windows.Automation.RangeValuePattern]::Pattern,[ref]$pattern)) {$s=$pattern.Current;$range=@{minimum=$s.Minimum;maximum=$s.Maximum;value=$s.Value;small_change=$s.SmallChange;large_change=$s.LargeChange;read_only=$s.IsReadOnly}}
   $pattern=$null
   if ($element.TryGetCurrentPattern([System.Windows.Automation.ScrollPattern]::Pattern,[ref]$pattern)) {$s=$pattern.Current;$scroll=@{vertical_percent=$s.VerticalScrollPercent;vertical_view_size=$s.VerticalViewSize;vertically_scrollable=$s.VerticallyScrollable;horizontal_percent=$s.HorizontalScrollPercent;horizontal_view_size=$s.HorizontalViewSize;horizontally_scrollable=$s.HorizontallyScrollable}}
   if ($c.Name.Length -gt 4096 -or ($null -ne $value -and $value.Length -gt 4096)) {throw 'Owned UI text exceeds bound'}
   $automation.Add([pscustomobject]@{name=$c.Name;value=$value;control=$c.ControlType.ProgrammaticName;automation_id=$c.AutomationId;class=$c.ClassName;enabled=$c.IsEnabled;offscreen=$c.IsOffscreen;rect=@($c.BoundingRectangle.X,$c.BoundingRectangle.Y,$c.BoundingRectangle.Width,$c.BoundingRectangle.Height);patterns=@($element.GetSupportedPatterns() | ForEach-Object {$_.ProgrammaticName});range=$range;scroll=$scroll})
  }
  if ($process.HasExited -or $viewRoot.Current.ProcessId -ne $process.Id) {throw 'View observation owner/lifetime changed'}
 }
 $automation.ToArray() | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output ($label+'-automation.json')) -Encoding UTF8
 $main=@($automation.ToArray() | Where-Object {$_.automation_id -eq 'MainWindow.main_widget_.main_tab_widget_'})
 $loading=@($automation.ToArray() | Where-Object {$_.automation_id -eq 'MainWindow.file_loading_widget_' -and !$_.offscreen})
 $record.loaded_profile_observed=($main.Count -eq 1 -and $main[0].enabled -and $loading.Count -eq 0)
 if (!$record.loaded_profile_observed) {throw 'Loaded profile readiness lost'}
 if (!@($rows | Where-Object {$_.print_window_succeeded -and $_.image_bytes -gt 0 -and $_.image_sha256}).Count) {throw 'Owned view image unavailable'}
 foreach ($row in $rows) {
  if ($row.print_window_succeeded -and $row.image_path) {
   $target=Join-Path $Output ($label+'-'+[IO.Path]::GetFileName($row.image_path))
   if (Test-Path -LiteralPath $target) {throw 'One-shot view image already exists'}
   Move-Item -LiteralPath $row.image_path -Destination $target;$row.image_path=$target
  }
 }
 $record.observations+=@([pscustomobject]@{label=$label;automation_path=(Join-Path $Output ($label+'-automation.json'));ui_elements=$automation.Count;window_observation=$rows})
 $record | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
}
'''
ps=replace(ps,'$process=$null;$watch=',functions+observe+'\n$process=$null;$watch=')
start=ps.index('    $record.navigation_history=@()')
end=ps.index('    $record.intentional_owned_stop=$true',start)
actions='''    $record.navigation_history=@();$record.observations=@()
    SelectOwnedTab 'EVENTS'
    SelectOwnedTab 'Event timing'
    ObserveOwnedView 'events'
    SelectOwnedTab 'Wavefront occupancy'
    ObserveOwnedView 'counter-graphs'
    GuardOwnedView
    $condition=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::AutomationIdProperty,'BUTTON_ID')
    $buttons=$root.FindAll([System.Windows.Automation.TreeScope]::Descendants,$condition)
    if ($buttons.Count -ne 1) {throw 'Observed Views button must remain unique'}
    $button=$buttons[0];$c=$button.Current
    if ($c.ProcessId -ne $process.Id -or !$c.IsEnabled -or $c.IsOffscreen -or $c.Name -ne 'Views' -or $c.ClassName -ne 'ArrowIconComboBox' -or $process.HasExited -or $root.Current.ProcessId -ne $process.Id) {throw 'Observed Views target owner/lifetime/state mismatch'}
    $r=$button.Current.BoundingRectangle
    $record.views_menu_request=@{hwnd=$root.Current.NativeWindowHandle;pid=$process.Id;rect=@($r.X,$r.Y,$r.Width,$r.Height)}
    $record | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
    $record.views_menu_message_result=[SiriusWindow]::ClickOwned($root.Current.NativeWindowHandle,[uint32]$process.Id,$r.X,$r.Y,$r.Width,$r.Height)
    Start-Sleep -Milliseconds 500
    GuardOwnedView
    ObserveOwnedView 'counter-view-menu'
    $record.navigation_completed=$true
    $record.window_observation=@($record.observations | ForEach-Object {$_.window_observation})
'''.replace('BUTTON_ID',button['automation_id'])
ps=ps[:start]+actions+ps[end:]
ps=replace(ps,"if (!$record.navigation_completed -or", "if ($record.observations.Count -ne 3 -or !$record.navigation_completed -or")
ps=ps.replace('Owned sealed-profile EVENTS and Wavefront occupancy tab selection and GUI observation only','Owned saved-profile actual Event timing, Wavefront counter graphs, Views menu and read-only pattern inventory')
save('profiler_counter_controls.ps1',ps)
b=(WORK/'run_profiler_counter_inventory.py').read_text().replace('profiler-counter-inventory','profiler-counter-controls').replace('profiler_counter_inventory.ps1','profiler_counter_controls.ps1').replace('counter-inventory-source-review.json','counter-controls-source-review.json')
b=replace(b,"paths += [trace,", "paths += [WORK/'prepare_counter_controls.py',WORK/'profiler-counter-inventory-bridge.json',WORK/'profiler-counter-inventory/automation.json',trace,")
b=replace(b,"['EVENTS','Wavefront occupancy']", "['EVENTS','Event timing','Wavefront occupancy']")
b=replace(b,"assert o['test_executable_sha256']", "assert len(o['observations'])==3 and [v['label'] for v in o['observations']]==['events','counter-graphs','counter-view-menu']\n  assert o['test_executable_sha256']")
save('run_profiler_counter_controls.py',b)
plan={'accepted_prior_inventory':seal(WORK/'profiler-counter-inventory-bridge.json'),'accepted_prior_automation':seal(WORK/'profiler-counter-inventory/automation.json'),'observed_views_button':button,'scope':'Current trace event inventory and actual existing graph/menu/pattern discovery; no event selection, graph region selection or counter attribution. All prior ownership and time/memory guards retained.'}
save('counter-controls-plan.json',json.dumps(plan,indent=2)+'\n')
print(json.dumps({'prepared':True,'scope':plan['scope']}))
