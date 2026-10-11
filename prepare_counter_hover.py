"""Use reversible real cursor hovers where synthetic messages returned no tooltips."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import ast,json,hashlib
W=Path(__file__).resolve().parent
def replace(s,o,n):assert s.count(o)==1,(o[:80],s.count(o));return s.replace(o,n)
def save(n,s):
 p=W/n;assert not p.exists()
 if n.endswith('.py'):ast.parse(s,filename=n)
 p.write_text(s)
prior=json.loads((W/'profiler-counter-event16-bridge.json').read_text());assert prior['accepted']
ps=(W/'profiler_counter_event.ps1').read_text().replace('profiler-counter-event16','profiler-counter-hover16').replace('profiler-counter-event22','profiler-counter-hover22')
ps=replace(ps,'public sealed class SiriusWindowInfo', (W/'owned_cursor.cs').read_text().replace('using System;\nusing System.Runtime.InteropServices;\n','')+'\npublic sealed class SiriusWindowInfo')
ps=replace(ps,'if (w<100 || y<100)', 'if (w<20 || y<20)')
ps=replace(ps,"    SelectOwnedTab 'Pipeline state'\n    ObserveOwnedView 'pipeline'\n",'')
old="     $record.counter_hovers+=@(@{id=$item.id;rect=@($r.X,$r.Y,$r.Width,$r.Height);message=[SiriusWindow]::PointerOwned($root.Current.NativeWindowHandle,[uint32]$process.Id,$r.X,$r.Y,$r.Width,$r.Height,$false)})"
new='''     if (!$record.Contains('cursor')) {
      $before=[SiriusCursor]::Read();$record.cursor=@{before=$before;last=$before;target=$before;foreground=[SiriusCursor]::Foreground();restored=$false;tool_pid=$process.Id;tool_start_utc=$record.test_start_time_utc}
     }
     $now=[SiriusCursor]::Read()
     if ($now[0] -ne $record.cursor.last[0] -or $now[1] -ne $record.cursor.last[1] -or [SiriusCursor]::Foreground() -ne $record.cursor.foreground) {throw 'External cursor/foreground change before hover'}
     $record.cursor.target=@([int][Math]::Floor($r.X+$r.Width/2),[int][Math]::Floor($r.Y+$r.Height/2))
     $record.cursor | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'cursor-state.json') -Encoding UTF8
     $record.cursor.last=[SiriusCursor]::HoverOwned($root.Current.NativeWindowHandle,[uint32]$process.Id,$r.X,$r.Y,$r.Width,$r.Height,$record.cursor.last[0],$record.cursor.last[1],$record.cursor.foreground)
     $record.cursor | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'cursor-state.json') -Encoding UTF8
     $record.counter_hovers+=@(@{id=$item.id;rect=@($r.X,$r.Y,$r.Width,$r.Height);cursor=$record.cursor.last;scope='Real cursor hover within owned client; original position restored unless externally changed; no click/input/activation'})'''
ps=replace(ps,old,new)
restore='''function RestoreOwnedCursor {
 if ($record.Contains('cursor') -and !$record.cursor.restored) {
  $c=$record.cursor
  $c.result=[SiriusCursor]::Restore($c.before[0],$c.before[1],$c.last[0],$c.last[1],$c.target[0],$c.target[1])
  $c.restored=$c.result -in @('restored','already-original','external-position-preserved')
  $c | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'cursor-state.json') -Encoding UTF8
  if (!$c.restored) {throw 'Cursor restoration unresolved'}
 }
}
'''
ps=replace(ps,'function FindOwnedId',restore+'\nfunction FindOwnedId')
ps=replace(ps,'    $record.intentional_owned_stop=$true', '    RestoreOwnedCursor\n    $record.intentional_owned_stop=$true')
ps=replace(ps,'finally {\n $record.owned_process_absent=', "finally {\n try {RestoreOwnedCursor} catch {$errors.Add('cursor-restoration: '+$_.Exception.Message)}\n $record.owned_process_absent=")
ps=replace(ps,'$record.observations.Count -ne 7','$record.observations.Count -ne 6')
save('profiler_counter_hover.ps1',ps)
b=(W/'run_profiler_counter_event.py').read_text().replace('profiler-counter-event16','profiler-counter-hover16').replace('profiler-counter-event22','profiler-counter-hover22').replace('profiler_counter_event.ps1','profiler_counter_hover.ps1').replace('counter-event-source-review.json','counter-hover-source-review.json')
b=replace(b,"event_name=event_plan['events'][action]", "event_name=event_plan['events'][action.replace('hover','event')]")
b=replace(b,"paths += [WORK/'counter-event-plan.json'", "paths += [WORK/'counter-hover-plan.json',WORK/'owned_cursor.cs',WORK/'restore_owned_cursor.ps1',WORK/'prepare_counter_hover.py',WORK/'profiler-counter-event16-bridge.json',WORK/'counter-event-plan.json'")
old='& "+quote(win(WORK/\'emergency_native_cleanup.ps1\'))+'
new='try {& "+quote(win(WORK/\'restore_owned_cursor.ps1\'))+" -Output "+quote(win(output))+" -CursorSource "+quote(win(WORK/\'owned_cursor.cs\'))+" -CursorHash "+quote(sha(WORK/\'owned_cursor.cs\'))+"} finally {& "+quote(win(WORK/\'emergency_native_cleanup.ps1\'))+'
b=replace(b,old,new)
b=replace(b,"+quote(win(WORK/(action+'-emergency.json')))\ncleanup=","+quote(win(WORK/(action+'-emergency.json')))+\"}\"\ncleanup=")
b=replace(b,"assert receipt['recorded_births_absent'] and not receipt['errors']", "assert receipt['recorded_births_absent'] and not receipt['errors']\n  assert read(output/'cursor-recovery.json')['pass']")
b=replace(b,"['EVENTS','Event timing','Pipeline state','Wavefront occupancy']", "['EVENTS','Event timing','Wavefront occupancy']")
b=replace(b,"len(o['observations'])==7", "len(o['observations'])==6")
b=replace(b,"['event','pipeline','region','counter-graphs','cache-tooltip','lds-tooltip','memory-tooltip']", "['event','region','counter-graphs','cache-tooltip','lds-tooltip','memory-tooltip']")
b=replace(b,"assert o['selected_event_confirmed']", "assert o['cursor']['restored'] and o['cursor']['result'] in ('restored','already-original','external-position-preserved')\n  assert o['selected_event_confirmed']")
save('run_profiler_counter_hover.py',b)
save('counter-hover-plan.json',json.dumps({'prior_synthetic_hover_accepted':True,'prior_synthetic_hover_numeric_tooltip_observed':False,'change':'Real cursor hover only when current and target positions already belong to owned visible client, no activation/clicks; exact original restoration with externally changed position preserved. Scoped outer recovery restores recorded token before tool cleanup.','scope':'Fresh saved-profile event/region/counter tooltip observations under unchanged60s/1GiB/120s guards; no acquisition, science or benefit.'},indent=2)+'\n')
print('Prepared bounded reversible owned-client cursor hovers.')
