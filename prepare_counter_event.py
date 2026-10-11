"""Observe two actually inventoried dispatches and explicit RGP counter regions."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import ast,hashlib,json
WORK=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def replace(s,o,n):
 assert s.count(o)==1,(o[:80],s.count(o));return s.replace(o,n)
def save(n,s):
 p=WORK/n;assert not p.exists()
 if n.endswith('.py'):ast.parse(s,filename=n)
 p.write_text(s)
prior=json.loads((WORK/'profiler-counter-controls-bridge.json').read_text());assert prior['accepted']
events=json.loads((WORK/'profiler-counter-controls/events-automation.json').read_text(encoding='utf-8-sig'))
graphs=json.loads((WORK/'profiler-counter-controls/counter-graphs-automation.json').read_text(encoding='utf-8-sig'))
names={'profiler-counter-event16':'16 vkCmdDispatch(3, 1, 1)','profiler-counter-event22':'22 vkCmdDispatch(6, 1, 1)'}
for name in names.values():
 matches=[v for v in events if v['name']==name and v['control']=='ControlType.TreeItem'];assert len(matches)==1 and matches[0]['enabled'] and not matches[0]['offscreen']
base=next(v['automation_id'].rsplit('.',1)[0] for v in graphs if v['automation_id'].endswith('.customize_views_'))
for name in ('wavefront_occupancy_view_','spm_view_Cache_','spm_view_LDS_','spm_view_Memory (%)_'):
 assert len([v for v in graphs if v['automation_id']==name])==1
container=base+'.customizable_widget_scroll_area_.qt_scrollarea_vcontainer'
assert len([v for v in graphs if v['automation_id']==container])==1
ps=(WORK/'profiler_counter_controls.ps1').read_text()
ps=replace(ps,"[ValidateSet('profiler-counter-controls')]", "[ValidateSet('profiler-counter-event16','profiler-counter-event22')]")
ps=replace(ps,'[string]$TraceHash)', '[string]$TraceHash,[Parameter(Mandatory=$true)][string]$EventName)')
pointer=r'''
 public static string PointerOwned(long hwnd, uint pid, double x, double y, double width, double height, bool drag) {
  if(double.IsNaN(x)||double.IsNaN(y)||double.IsNaN(width)||double.IsNaN(height)||double.IsInfinity(x)||double.IsInfinity(y)||double.IsInfinity(width)||double.IsInfinity(height)||width<20||height<20||width>4096||height>4096||x<int.MinValue||y<int.MinValue||x+width>int.MaxValue||y+height>int.MaxValue) throw new Exception("Observed graph rectangle is invalid");
  IntPtr h=new IntPtr(hwnd);uint actual;RECT client;
  if(GetWindowThreadProcessId(h,out actual)==0||actual!=pid||!IsWindowVisible(h)||!GetClientRect(h,out client)||client.R-client.L>4096||client.B-client.T>4096) throw new Exception("Graph HWND/client not visibly owned");
  POINT first=new POINT{X=(int)Math.Floor(x),Y=(int)Math.Floor(y)},last=new POINT{X=(int)Math.Ceiling(x+width)-1,Y=(int)Math.Ceiling(y+height)-1};
  if(!ScreenToClient(h,ref first)||!ScreenToClient(h,ref last)||first.X<client.L||first.Y<client.T||last.X>=client.R||last.Y>=client.B) throw new Exception("Graph is clipped outside owned client");
  int left=first.X+(last.X-first.X)/10,right=first.X+9*(last.X-first.X)/10,cy=(first.Y+last.Y)/2;
  if(left<0||cy<0||right>4095||cy>4095) throw new Exception("Graph client coordinates exceed bound");
  Action<uint,uint,int> send=(message,buttons,cx)=>{UIntPtr result;if(GetWindowThreadProcessId(h,out actual)==0||actual!=pid) throw new Exception("Graph message owner changed");if(SendMessageTimeout(h,message,new UIntPtr(buttons),new IntPtr((cy<<16)|cx),0x22,1000,out result)==IntPtr.Zero) throw new Exception("Graph message failed or timed out");if(GetWindowThreadProcessId(h,out actual)==0||actual!=pid) throw new Exception("Graph message owner changed after send");};
  if(drag) {send(0x200,0,left);send(0x201,1,left);for(int i=1;i<=8;i++) send(0x200,1,left+(right-left)*i/8);send(0x202,0,right);}
  else {send(0x200,0,(first.X+last.X)/2);}
  return (drag?"Owned graph drag":"Owned graph hover")+" at client x="+left+".."+right+",y="+cy+"; no global cursor/input API";
 }
'''
ps=replace(ps,' public static SiriusWindowInfo[] Capture',pointer+' public static SiriusWindowInfo[] Capture')
functions=r'''
function FindOwnedId([string]$id) {
 GuardOwnedView
 if ($root.Current.ProcessId -ne $process.Id) {throw 'Owned root identity changed'}
 $condition=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::AutomationIdProperty,$id)
 $items=$root.FindAll([System.Windows.Automation.TreeScope]::Descendants,$condition)
 if ($items.Count -ne 1) {throw ('Expected unique owned id: '+$id)}
 $item=$items[0];$c=$item.Current
 if ($c.ProcessId -ne $process.Id -or !$c.IsEnabled -or $c.IsOffscreen) {throw 'Owned id state mismatch'}
 return $item
}
function ClickOwnedElement($item) {
 GuardOwnedView
 $c=$item.Current;$r=$c.BoundingRectangle
 if ($c.ProcessId -ne $process.Id -or !$c.IsEnabled -or $c.IsOffscreen -or $root.Current.ProcessId -ne $process.Id) {throw 'Owned click target state mismatch'}
 $record.last_owned_click=@{id=$c.AutomationId;name=$c.Name;rect=@($r.X,$r.Y,$r.Width,$r.Height)}
 $record | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 $record.last_owned_click_result=[SiriusWindow]::ClickOwned($root.Current.NativeWindowHandle,[uint32]$process.Id,$r.X,$r.Y,$r.Width,$r.Height)
 Start-Sleep -Milliseconds 500
 GuardOwnedView
}
'''
ps=replace(ps,'function ObserveOwnedView',functions+'\nfunction ObserveOwnedView')
start=ps.index('    $record.navigation_history=@();$record.observations=@()');end=ps.index('    $record.intentional_owned_stop=$true',start)
actions=r'''    $record.navigation_history=@();$record.observations=@()
    $record.requested_event=$EventName
    SelectOwnedTab 'EVENTS'
    SelectOwnedTab 'Event timing'
    $nameCondition=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::NameProperty,$EventName)
    $typeCondition=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::ControlTypeProperty,[System.Windows.Automation.ControlType]::TreeItem)
    $eventCondition=New-Object System.Windows.Automation.AndCondition($nameCondition,$typeCondition)
    $matches=$root.FindAll([System.Windows.Automation.TreeScope]::Descendants,$eventCondition)
    if ($matches.Count -ne 1) {throw 'Current observed dispatch is not unique'}
    ClickOwnedElement $matches[0]
    $matches=$root.FindAll([System.Windows.Automation.TreeScope]::Descendants,$eventCondition);$pattern=$null
    if ($matches.Count -ne 1 -or $matches[0].Current.ProcessId -ne $process.Id -or !$matches[0].TryGetCurrentPattern([System.Windows.Automation.SelectionItemPattern]::Pattern,[ref]$pattern) -or !$pattern.Current.IsSelected) {throw 'Fresh owned dispatch selection missing'}
    $record.selected_event_confirmed=$true
    ObserveOwnedView 'event'
    SelectOwnedTab 'Pipeline state'
    ObserveOwnedView 'pipeline'
    SelectOwnedTab 'Wavefront occupancy'
    $selected=FindOwnedId 'BASE.selected_event_name_'
    $record.wavefront_selected_name=$selected.Current.Name
    if ($record.wavefront_selected_name -ne $EventName) {throw 'Wavefront selected event identity mismatch'}
    $zoom=FindOwnedId 'BASE.zoom_to_selection_'
    ClickOwnedElement $zoom
    $graph=FindOwnedId 'wavefront_occupancy_view_';$r=$graph.Current.BoundingRectangle
    $record.region_request=@{event=$EventName;graph_id='wavefront_occupancy_view_';rect=@($r.X,$r.Y,$r.Width,$r.Height);requested_horizontal_fraction=@(0.1,0.9);scope='Explicit interior of event-zoomed visible ruler; actual returned region endpoints/duration authoritative. No exclusive shader attribution.'}
    $record | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
    $record.region_message_result=[SiriusWindow]::PointerOwned($root.Current.NativeWindowHandle,[uint32]$process.Id,$r.X,$r.Y,$r.Width,$r.Height,$true)
    Start-Sleep -Milliseconds 500
    GuardOwnedView
    ObserveOwnedView 'region'
    $scrollContainer=FindOwnedId 'CONTAINER'
    $scrollType=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::ControlTypeProperty,[System.Windows.Automation.ControlType]::ScrollBar)
    $bars=$scrollContainer.FindAll([System.Windows.Automation.TreeScope]::Descendants,$scrollType);$range=$null
    if ($bars.Count -ne 1 -or $bars[0].Current.ProcessId -ne $process.Id -or !$bars[0].TryGetCurrentPattern([System.Windows.Automation.RangeValuePattern]::Pattern,[ref]$range) -or $range.Current.IsReadOnly) {throw 'Observed graph vertical scrollbar unavailable'}
    $max=$range.Current.Maximum
    if ([double]::IsNaN($max) -or [double]::IsInfinity($max) -or $max -lt 0 -or $max -gt 4096) {throw 'Counter vertical scroll bound invalid'}
    $record.counter_scroll_request=@{container='CONTAINER';before=$range.Current.Value;maximum=$max}
    $range.SetValue($max);Start-Sleep -Milliseconds 300
    GuardOwnedView
    if ($bars[0].Current.ProcessId -ne $process.Id -or $range.Current.Value -ne $max) {throw 'Owned counter scrollbar did not reach maximum'}
    ObserveOwnedView 'counter-graphs'
    $record.counter_hovers=@()
    foreach ($item in @(@{id='spm_view_Cache_';label='cache-tooltip'},@{id='spm_view_LDS_';label='lds-tooltip'},@{id='spm_view_Memory (%)_';label='memory-tooltip'})) {
     GuardOwnedView
     $graph=FindOwnedId $item.id;$r=$graph.Current.BoundingRectangle
     $record.counter_hovers+=@(@{id=$item.id;rect=@($r.X,$r.Y,$r.Width,$r.Height);message=[SiriusWindow]::PointerOwned($root.Current.NativeWindowHandle,[uint32]$process.Id,$r.X,$r.Y,$r.Width,$r.Height,$false)})
     Start-Sleep -Milliseconds 1100
     GuardOwnedView
     ObserveOwnedView $item.label
    }
    $record.navigation_completed=$true
    $record.window_observation=@($record.observations | ForEach-Object {$_.window_observation})
'''.replace('BASE',base).replace('CONTAINER',container)
ps=ps[:start]+actions+ps[end:]
ps=replace(ps,'$record.observations.Count -ne 3','$record.observations.Count -ne 7')
ps=ps.replace('Owned saved-profile actual Event timing, Wavefront counter graphs, Views menu and read-only pattern inventory','Owned current-profile selected dispatch/pipeline, explicit graph-region and counter hover observations; actual returned data needed for semantics')
save('profiler_counter_event.ps1',ps)
b=(WORK/'run_profiler_counter_controls.py').read_text().replace("('profiler-counter-controls',)","('profiler-counter-event16','profiler-counter-event22')").replace('profiler_counter_controls.ps1','profiler_counter_event.ps1').replace('counter-controls-source-review.json','counter-event-source-review.json')
b=replace(b,"control_plan=read(WORK/'counter-controls-plan.json')", "event_plan=read(WORK/'counter-event-plan.json');event_name=event_plan['events'][action]\ncontrol_plan=read(WORK/'counter-controls-plan.json')")
b=replace(b,"paths += [WORK/'counter-controls-plan.json'", "paths += [WORK/'counter-event-plan.json',WORK/'prepare_counter_event.py',WORK/'profiler-counter-controls-bridge.json',WORK/'profiler-counter-controls/events-automation.json',WORK/'profiler-counter-controls/counter-graphs-automation.json',WORK/'counter-controls-plan.json'")
b=replace(b,"'-TraceHash',sha(trace)]", "'-TraceHash',sha(trace),'-EventName',event_name]")
b=replace(b,"['EVENTS','Event timing','Wavefront occupancy']", "['EVENTS','Event timing','Pipeline state','Wavefront occupancy']")
b=replace(b,"assert len(o['observations'])==3 and [v['label'] for v in o['observations']]==['events','counter-graphs','counter-view-menu']", "assert o['selected_event_confirmed'] and o['requested_event']==event_name and o['wavefront_selected_name']==event_name\n  assert len(o['observations'])==7 and [v['label'] for v in o['observations']]==['event','pipeline','region','counter-graphs','cache-tooltip','lds-tooltip','memory-tooltip']")
b=replace(b,"event_plan=read(WORK/'counter-event-plan.json');event_name=event_plan['events'][action]", "event_plan=read(WORK/'counter-event-plan.json');event_name=event_plan['events'][action]\nfor name,v in event_plan['prior_seals'].items():\n p=WORK/name;assert {'bytes':p.stat().st_size,'sha256':sha(p)}==v\ncurrent=read(WORK/'profiler-counter-controls-bridge.json');assert current['accepted'] and current['source_revision']==head")
save('run_profiler_counter_event.py',b)
plan={'events':names,'prior_seals':{n:seal(WORK/n) for n in ('profiler-counter-controls-bridge.json','profiler-counter-controls/events-automation.json','profiler-counter-controls/counter-graphs-automation.json')},'region_semantics':'Documented explicit WFO graph drag after selected-event zoom; requested fractions are geometric, actual returned endpoints/duration/tooltip required. Zoom alone is not aggregate. Synthetic owned mouse messages may fail to produce a tooltip; absence must be reported.','scope':'One owned event per action under original viewer/bridge bounds; no scientific execution, raw timestamp overlay, source binding or performance benefit assumed.'}
save('counter-event-plan.json',json.dumps(plan,indent=2)+'\n')
print(json.dumps({'prepared':True,'events':names}))
