"""Send bounded client mouse messages only to the owned RGP HWND; inspect effect."""
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
ps = (WORK / 'profiler_trace_invoke.ps1').read_text().replace(
    'profiler-trace-invoke', 'profiler-trace-owned-click')
native = """ [StructLayout(LayoutKind.Sequential)] private struct POINT {public int X,Y;}
 [DllImport("user32.dll", SetLastError=true)] private static extern bool ScreenToClient(IntPtr h, ref POINT p);
 [DllImport("user32.dll", SetLastError=true)] private static extern bool GetClientRect(IntPtr h, out RECT r);
 [DllImport("user32.dll", CharSet=CharSet.Unicode, SetLastError=true)] private static extern IntPtr SendMessageTimeout(IntPtr h, uint m, UIntPtr w, IntPtr l, uint f, uint t, out UIntPtr result);
 public static string ClickOwned(long hwnd, uint pid, double x, double y, double width, double height) {
  if(double.IsNaN(x)||double.IsNaN(y)||double.IsNaN(width)||double.IsNaN(height)||double.IsInfinity(x)||double.IsInfinity(y)||double.IsInfinity(width)||double.IsInfinity(height)||width<=0||height<=0||width>4096||height>4096||x<int.MinValue||y<int.MinValue||x+width>int.MaxValue||y+height>int.MaxValue) throw new Exception("Observed target rectangle is invalid");
  IntPtr h=new IntPtr(hwnd);uint actual;
  if(GetWindowThreadProcessId(h,out actual)==0||actual!=pid||!IsWindowVisible(h)) throw new Exception("Mouse-message HWND is not visibly owned");
  RECT client;if(!GetClientRect(h,out client)||client.R-client.L>4096||client.B-client.T>4096) throw new Exception("Owned client rectangle unavailable or excessive");
  POINT first=new POINT{X=(int)Math.Floor(x),Y=(int)Math.Floor(y)},last=new POINT{X=(int)Math.Ceiling(x+width)-1,Y=(int)Math.Ceiling(y+height)-1};
  if(!ScreenToClient(h,ref first)||!ScreenToClient(h,ref last)||first.X<client.L||first.Y<client.T||last.X>=client.R||last.Y>=client.B) throw new Exception("Observed target is outside owned client rectangle");
  int cx=(first.X+last.X)/2,cy=(first.Y+last.Y)/2;
  if(cx<0||cy<0||cx>4095||cy>4095) throw new Exception("Owned mouse-message coordinates exceed bound");
  IntPtr coords=new IntPtr((cy<<16)|cx);UIntPtr result;
  if(GetWindowThreadProcessId(h,out actual)==0||actual!=pid) throw new Exception("Mouse-message owner changed before down");
  if(SendMessageTimeout(h,0x201,new UIntPtr(1),coords,0x22,1000,out result)==IntPtr.Zero) throw new Exception("Owned button-down message failed or timed out");
  if(GetWindowThreadProcessId(h,out actual)==0||actual!=pid) throw new Exception("Mouse-message owner changed before up");
  if(SendMessageTimeout(h,0x202,UIntPtr.Zero,coords,0x22,1000,out result)==IntPtr.Zero) throw new Exception("Owned button-up message failed or timed out");
  if(GetWindowThreadProcessId(h,out actual)==0||actual!=pid) throw new Exception("Mouse-message owner changed after up");
  return "Owned HWND WM_LBUTTONDOWN/UP completed at client "+cx+","+cy+"; no global input/cursor/focus API used";
 }
"""
ps = replace(ps, ' public static SiriusWindowInfo[] Capture(uint pid, string output) {', native+' public static SiriusWindowInfo[] Capture(uint pid, string output) {')
ps = replace(ps, "     $selection=$null\n     if (!$event.TryGetCurrentPattern([System.Windows.Automation.InvokePattern]::Pattern,[ref]$selection)) {throw 'Dispatch selection pattern unavailable'}\n", '')
ps = replace(ps, '     $selection.Invoke();Start-Sleep -Milliseconds 500', "     $r=$event.Current.BoundingRectangle\n     $record.owned_message_target=[ordered]@{hwnd=$rows[0].handle;pid=$process.Id;event_name=$eventName;rect=@($r.X,$r.Y,$r.Width,$r.Height)}\n     $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8\n     if ($process.HasExited -or $event.Current.ProcessId -ne $process.Id -or $root.Current.ProcessId -ne $process.Id) {throw 'Owned mouse-message target lifetime changed'}\n     $record.owned_mouse_message_result=[SiriusWindow]::ClickOwned($rows[0].handle,[uint32]$process.Id,$r.X,$r.Y,$r.Width,$r.Height)\n     Start-Sleep -Milliseconds 500")
ps = replace(ps, "requested_action='InvokePattern.Invoke'", "requested_action='Owned HWND client WM_LBUTTONDOWN/UP'")
ps = replace(ps, 'Owned saved-profile tabs and one observed dispatch InvokePattern request/result observation only; its effect is not assumed;', 'Owned saved-profile tabs and one bounded HWND client mouse-message request/result observation only; no global input API or assumed effect;')
save('profiler_trace_owned_click.ps1', ps)
b = (WORK / 'run_profiler_trace_invoke.py').read_text().replace(
    'profiler-trace-invoke', 'profiler-trace-owned-click').replace(
    'profiler_trace_invoke.ps1', 'profiler_trace_owned_click.ps1').replace(
    'trace-invoke-source-review.json', 'trace-owned-click-source-review.json').replace(
    'prepare_trace_invoke.py', 'prepare_trace_owned_click.py')
b = replace(b, 'Owned saved-profile tabs and one dispatch InvokePattern request/result observation only; no effect assumed.', 'Owned saved-profile tabs and one bounded owned-HWND client mouse-message request/result observation only; no global input or effect assumed.')
b = replace(b, "WORK/'profiler-trace-selection/owner.json'", "WORK/'profiler-trace-selection/owner.json',WORK/'profiler-trace-invoke-bridge.json',WORK/'profiler-trace-invoke/owner.json'")
save('run_profiler_trace_owned_click.py', b)
