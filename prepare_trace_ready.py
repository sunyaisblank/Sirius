"""Wait for the existing profile loader to finish; no new capture or science."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
WORK=Path(__file__).resolve().parent
def replace(s,old,new):
 assert s.count(old)==1,(old,s.count(old));return s.replace(old,new)
def save(n,s):
 p=WORK/n;assert not p.exists();p.write_text(s)
ps=(WORK/'profiler_trace_ui.ps1').read_text()
ps=replace(ps,"'profiler-trace-ui'","'profiler-trace-ready'")
ps=replace(ps,'polling_guard_seconds=30;blocked_call_bridge_guard_seconds=90','polling_guard_seconds=60;blocked_call_bridge_guard_seconds=120')
ps=replace(ps,'$watch.Elapsed.TotalSeconds -gt 30','$watch.Elapsed.TotalSeconds -gt 60')
ps=replace(ps,'    $record.owned_ui_elements=$automation.Count', """    $record.owned_ui_elements=$automation.Count
    $main=@($automation.ToArray() | Where-Object {$_.automation_id -eq 'MainWindow.main_widget_.main_tab_widget_'})
    $loading=@($automation.ToArray() | Where-Object {$_.automation_id -eq 'MainWindow.file_loading_widget_' -and !$_.offscreen})
    $record.loaded_profile_observed=($main.Count -eq 1 -and $main[0].enabled -and $loading.Count -eq 0)
    if (!$record.loaded_profile_observed) {Start-Sleep -Milliseconds 200;continue}
    # The first screenshot may precede the loader transition; recapture now.
    $rows=@([SiriusWindow]::Capture([uint32]$process.Id,$Output))
    if ($process.HasExited) {throw 'Owned profiler exited during loaded-profile observation'}
    $record.window_observation=$rows""")
ps=replace(ps,"if ($record.status -ne 'observed-and-stopped'", "if (!$record.loaded_profile_observed -or $record.status -ne 'observed-and-stopped'")
save('profiler_trace_ready.ps1',ps)
b=(WORK/'run_profiler_trace_ui.py').read_text()
b=b.replace('profiler-trace-ui','profiler-trace-ready').replace('profiler_trace_ui.ps1','profiler_trace_ready.ps1').replace('trace-ui-source-review.json','trace-ready-source-review.json')
b=replace(b,"WORK/'prepare_trace_ui.py'", "WORK/'prepare_trace_ready.py',WORK/'profiler-trace-ui-bridge.json'")
# The unchanged exact-birth emergency helper still has its own 60s bound.
b=replace(b,"'deadline_seconds':90", "'deadline_seconds':120")
b=replace(b,'remaining=90-(time.monotonic()-start)','remaining=120-(time.monotonic()-start)')
b=replace(b,'raise subprocess.TimeoutExpired(args,90)','raise subprocess.TimeoutExpired(args,120)')
b=replace(b,"o['status']=='observed-and-stopped'", "o['loaded_profile_observed'] and o['status']=='observed-and-stopped'")
save('run_profiler_trace_ready.py',b)
