"""Derive an owned-window inspection from the reviewed one-shot help owner."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib
work = Path(__file__).resolve().parent
ps = (work / 'tool_help.ps1').read_text()
bridge = (work / 'run_tool_help.py').read_text()
assert hashlib.sha256(ps.encode()).hexdigest() == 'b5f166858f531c7a0edf565ad7cd22385078f5f6f3c95f75501e2f64c3bdc5cd'
assert hashlib.sha256(bridge.encode()).hexdigest() == 'a63df2d251fa4c7f16864daa7d9f43e3c498807acc9ca9dfa75292d35b521c05'
def change(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new)
ps = change(ps, "[ValidateSet('panel-help','panel-version','profiler-help')]", "[ValidateSet('profiler-ui')]")
ps = change(ps, "$tool=if ($Action -eq 'profiler-help') {'RadeonGPUProfiler.exe'} else {'RadeonDeveloperPanelCLI.exe'}\n$argument=if ($Action -eq 'panel-version') {'--version'} else {'--help'}", "$tool='RadeonGPUProfiler.exe'")
ps = change(ps, 'Add-Type -TypeDefinition', 'Add-Type -ReferencedAssemblies System.Drawing -TypeDefinition')
ps = change(ps, '\nusing System.Runtime.InteropServices;', '\nusing System.Runtime.InteropServices;\nusing System.Collections.Generic;\nusing System.Drawing;\nusing System.Drawing.Imaging;\nusing System.IO;\nusing System.Text;')
window_code = r'''
public sealed class SiriusWindowInfo {
 public long handle; public uint pid; public string title; public int width; public int height;
 public bool print_window_succeeded; public string image_path;
}
public static class SiriusWindow {
 [StructLayout(LayoutKind.Sequential)] private struct RECT {public int L,T,R,B;}
 private delegate bool Callback(IntPtr h, IntPtr data);
 [DllImport("user32.dll")] private static extern bool EnumWindows(Callback cb, IntPtr data);
 [DllImport("user32.dll")] private static extern uint GetWindowThreadProcessId(IntPtr h, out uint pid);
 [DllImport("user32.dll")] private static extern bool IsWindowVisible(IntPtr h);
 [DllImport("user32.dll")] private static extern bool GetWindowRect(IntPtr h, out RECT r);
 [DllImport("user32.dll", CharSet=CharSet.Unicode)] private static extern int GetWindowText(IntPtr h, StringBuilder s, int count);
 [DllImport("user32.dll")] private static extern bool PrintWindow(IntPtr h, IntPtr dc, uint flags);
 public static SiriusWindowInfo[] Capture(uint pid, string output) {
  var rows = new List<SiriusWindowInfo>();
  bool enumerated = EnumWindows(delegate(IntPtr h, IntPtr data) {
   uint actual; GetWindowThreadProcessId(h, out actual);
   if (actual != pid || !IsWindowVisible(h)) return true;
   RECT r; if (!GetWindowRect(h, out r)) throw new Exception("Owned window rectangle unavailable");
   int w=r.R-r.L, y=r.B-r.T;
   if (w<100 || y<100) return true;
   if (w>4096 || y>4096) throw new Exception("Owned window image exceeds diagnostic dimensions");
   var s=new StringBuilder(1024); GetWindowText(h, s, s.Capacity);
   var row=new SiriusWindowInfo {handle=h.ToInt64(),pid=actual,title=s.ToString(),width=w,height=y};
   using(var bitmap=new Bitmap(w,y)) {
    using(var graphics=Graphics.FromImage(bitmap)) {
     IntPtr dc=graphics.GetHdc();
     try {row.print_window_succeeded=PrintWindow(h,dc,2);} finally {graphics.ReleaseHdc(dc);}
    }
    if(row.print_window_succeeded) {
     row.image_path=Path.Combine(output,"owned-window-"+h.ToInt64().ToString("x")+".png");
     bitmap.Save(row.image_path,ImageFormat.Png);
    }
   }
   rows.Add(row); return true;
  },IntPtr.Zero);
  if(!enumerated) throw new Exception("Owned window enumeration failed");
  return rows.ToArray();
 }
}
'''
ps = change(ps, "\n'@\n$existing=", window_code + "\n'@\n$existing=")
ps = change(ps, 'arguments=@($argument)', 'arguments=@();intentional_owned_stop=$false;window_observation=@()')
ps = change(ps, ' -ArgumentList $argument', '')
ps = change(ps, '  Start-Sleep -Milliseconds 50', '''  if ($watch.Elapsed.TotalSeconds -ge 3) {
   $rows=@([SiriusWindow]::Capture([uint32]$process.Id,$Output))
   if ($rows.Count) {
    if ($process.HasExited) {throw 'Owned profiler exited during window observation'}
    $record.window_observation=$rows;$record.intentional_owned_stop=$true
    $process.Kill();if (!$process.WaitForExit(5000)) {throw 'Observed owned profiler did not terminate'};break
   }
  }
  Start-Sleep -Milliseconds 50''')
ps = change(ps, "if ($record.status -ne 'failed') {$record.status=if ($record.outer_stop) {'invalidated'} else {'completed'}}", "if ($record.status -ne 'failed') {$record.status=if ($record.outer_stop) {'invalidated'} elseif ($record.intentional_owned_stop -and $record.window_observation.Count) {'observed-and-stopped'} else {'no-window-observation'}}")
ps = change(ps, "if ($record.status -ne 'completed' -or $record.returncode -ne 0) {exit 1}", "if ($record.status -ne 'observed-and-stopped' -or !$record.intentional_owned_stop -or !$record.window_observation.Count) {exit 1}")
ps = ps.replace('Pinned official executable help/version only; no service/capture/scientific result or developer-mode/clock/driver experiment authorized.', 'Owned profiler GUI feasibility observation only. Deliberate owned stop after its visible windows are recorded; no service/capture/scientific result or clock/driver changes.')
ps = ps.replace('Early exact tool-help observer birth; no service/capture/scientific launch authorized.', 'Early exact GUI observer birth; no service/capture/scientific launch authorized.')
bridge = change(bridge, "('panel-help','panel-version','profiler-help')", "('profiler-ui',)")
bridge = change(bridge, "tool='RadeonGPUProfiler.exe' if action=='profiler-help' else 'RadeonDeveloperPanelCLI.exe'", "tool='RadeonGPUProfiler.exe'")
bridge = bridge.replace('tool-help-source-review.json', 'gui-inspection-source-review.json').replace('run_tool_help.py', 'run_profiler_ui.py').replace('tool_help.ps1', 'profiler_ui.ps1')
bridge = change(bridge, "paths += [ROOT/v['path'] for v in extraction['files']]", "gui=read(WORK/'tools-gui-extraction.json');paths += [WORK/'tools-gui-extraction.json']\npaths += [ROOT/v['path'] for v in extraction['files']+gui['files']]")
bridge = change(bridge, "for v in extraction['files']:assert seals", "for v in extraction['files']+gui['files']:assert seals")
bridge = change(bridge, "o['status']=='completed' and o['returncode']==0", "o['status']=='observed-and-stopped' and o['intentional_owned_stop'] and o['window_observation']")
bridge = bridge.replace('Tool help/version only. No service/science/device capture;', 'Owned GUI window observation only. Intentional owned stop after observation; no service/science/device capture;')
for name, data in [('profiler_ui.ps1', ps), ('run_profiler_ui.py', bridge)]:
    with (work / name).open('x') as f: f.write(data)
    print(name, len(data.encode()), hashlib.sha256(data.encode()).hexdigest())
