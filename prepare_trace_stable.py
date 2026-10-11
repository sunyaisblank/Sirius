"""Poll stable loader controls before reading the changing profile subtree."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
WORK=Path(__file__).resolve().parent
def replace(s,old,new):
 assert s.count(old)==1,(old,s.count(old));return s.replace(old,new)
def save(n,s):
 p=WORK/n;assert not p.exists();p.write_text(s)
ps=(WORK/'profiler_trace_ready.ps1').read_text().replace("'profiler-trace-ready'","'profiler-trace-stable'")
ps=replace(ps,' $nativeHandle=$process.Handle'," $nativeHandle=$process.Handle\n $record.ui_unavailable_pid_zero=0")
ps=replace(ps,"   if ($rows.Count) {", """   if ($rows.Count) {
    $readyRoots=0
    foreach ($row in $rows) {
     if ($process.HasExited) {throw 'Owned profiler exited before loader poll'}
     $root=[System.Windows.Automation.AutomationElement]::FromHandle([IntPtr]$row.handle)
     if ($root.Current.ProcessId -ne $process.Id) {throw 'Loader root owner mismatch'}
     $mainCondition=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::AutomationIdProperty,'MainWindow.main_widget_.main_tab_widget_')
     $main=$root.FindFirst([System.Windows.Automation.TreeScope]::Descendants,$mainCondition)
     $loadingCondition=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::AutomationIdProperty,'MainWindow.file_loading_widget_')
     $loading=$root.FindFirst([System.Windows.Automation.TreeScope]::Descendants,$loadingCondition)
     if ($main -and $main.Current.ProcessId -eq $process.Id -and $main.Current.IsEnabled -and (!$loading -or ($loading.Current.ProcessId -eq $process.Id -and $loading.Current.IsOffscreen))) {$readyRoots++}
     if ($process.HasExited -or $root.Current.ProcessId -ne $process.Id) {throw 'Loader poll owner/lifetime changed'}
    }
    if ($readyRoots -ne 1) {Start-Sleep -Milliseconds 200;continue}
""")
ps=replace(ps,"      if ($c.ProcessId -ne $process.Id -or $process.HasExited) {throw 'UI child owner/lifetime mismatch'}", "      if ($process.HasExited) {throw 'Owned profiler exited during ready UI observation'}\n      if ($c.ProcessId -eq 0) {$record.ui_unavailable_pid_zero++;continue}\n      if ($c.ProcessId -ne $process.Id) {throw 'UI child foreign owner mismatch'}")
ps=replace(ps,' $record.owned_process_absent=($null -eq $process)', ' $record.owned_process_absent=($null -eq $process)\n $record.tool_already_exited_on_finally_entry=$null;$record.exception_cleanup_kill_requested=$false')
ps=replace(ps,"  try { if (!$process.HasExited) {$process.Kill();if (!$process.WaitForExit(5000)) {throw 'Owned help tool did not terminate'}} }", "  try {$record.tool_already_exited_on_finally_entry=$process.HasExited;if (!$record.tool_already_exited_on_finally_entry) {$record.exception_cleanup_kill_requested=$true;$process.Kill();if (!$process.WaitForExit(5000)) {throw 'Owned help tool did not terminate'}} }")
ps=replace(ps,"  try {$record.owned_process_absent=$process.HasExited}", "  try {$record.owned_process_absent=$process.HasExited;if ($record.owned_process_absent) {$record.terminal_returncode=[SiriusOwnedExit]::Read($nativeHandle);$record.terminal_exit_code_source='Held native handle after terminal wait, including failed observation'}}")
save('profiler_trace_stable.ps1',ps)
b=(WORK/'run_profiler_trace_ready.py').read_text().replace('profiler-trace-ready','profiler-trace-stable').replace('profiler_trace_ready.ps1','profiler_trace_stable.ps1').replace('trace-ready-source-review.json','trace-stable-source-review.json').replace('prepare_trace_ready.py','prepare_trace_stable.py')
b=replace(b,"o['loaded_profile_observed'] and", "o['terminal_returncode']==2**32-1 and o['loaded_profile_observed'] and")
save('run_profiler_trace_stable.py',b)
