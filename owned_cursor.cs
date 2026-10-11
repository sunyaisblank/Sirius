using System;
using System.Runtime.InteropServices;
public static class SiriusCursor {
 [StructLayout(LayoutKind.Sequential)] private struct POINT {public int X,Y;}
 [StructLayout(LayoutKind.Sequential)] private struct RECT {public int L,T,R,B;}
 [DllImport("user32.dll",SetLastError=true)] private static extern bool GetCursorPos(out POINT p);
 [DllImport("user32.dll",SetLastError=true)] private static extern bool SetCursorPos(int x,int y);
 [DllImport("user32.dll")] private static extern IntPtr WindowFromPoint(POINT p);
 [DllImport("user32.dll")] private static extern IntPtr GetAncestor(IntPtr h,uint flags);
 [DllImport("user32.dll")] private static extern uint GetWindowThreadProcessId(IntPtr h,out uint pid);
 [DllImport("user32.dll")] private static extern bool IsWindowVisible(IntPtr h);
 [DllImport("user32.dll")] private static extern bool ScreenToClient(IntPtr h,ref POINT p);
 [DllImport("user32.dll")] private static extern bool GetClientRect(IntPtr h,out RECT r);
 [DllImport("user32.dll")] private static extern IntPtr GetForegroundWindow();
 public static int[] Read() {POINT p;if(!GetCursorPos(out p)) throw new Exception("Cursor position unavailable");return new int[]{p.X,p.Y};}
 public static long Foreground() {return GetForegroundWindow().ToInt64();}
 private static void OwnedPoint(IntPtr h,uint pid,POINT p) {
  uint actual;RECT r;POINT local=p;IntPtr target=WindowFromPoint(p);
  if(GetWindowThreadProcessId(h,out actual)==0||actual!=pid||!IsWindowVisible(h)||target==IntPtr.Zero||GetWindowThreadProcessId(target,out actual)==0||actual!=pid||GetAncestor(target,2)!=h||!ScreenToClient(h,ref local)||!GetClientRect(h,out r)||local.X<r.L||local.Y<r.T||local.X>=r.R||local.Y>=r.B) throw new Exception("Cursor point is not in the owned visible client");
 }
 public static int[] HoverOwned(long hwnd,uint pid,double x,double y,double w,double height,int expectedX,int expectedY,long expectedForeground) {
  if(double.IsNaN(x)||double.IsNaN(y)||double.IsNaN(w)||double.IsNaN(height)||double.IsInfinity(x)||double.IsInfinity(y)||double.IsInfinity(w)||double.IsInfinity(height)||w<20||height<20||w>4096||height>4096||x<int.MinValue||y<int.MinValue||x+w>int.MaxValue||y+height>int.MaxValue) throw new Exception("Cursor graph rectangle invalid");
  IntPtr h=new IntPtr(hwnd);POINT before;if(!GetCursorPos(out before)) throw new Exception("Cursor unavailable before hover");
  // Respect the shared-cursor contract: its current and target positions must
  // already belong to this exact visible application client. No activation.
  OwnedPoint(h,pid,before);
  POINT first=new POINT{X=(int)Math.Floor(x),Y=(int)Math.Floor(y)},last=new POINT{X=(int)Math.Ceiling(x+w)-1,Y=(int)Math.Ceiling(y+height)-1};
  OwnedPoint(h,pid,first);OwnedPoint(h,pid,last);
  POINT target=new POINT{X=(int)Math.Floor(x+w/2),Y=(int)Math.Floor(y+height/2)};OwnedPoint(h,pid,target);
  POINT current;if(!GetCursorPos(out current)||current.X!=expectedX||current.Y!=expectedY||GetForegroundWindow().ToInt64()!=expectedForeground) throw new Exception("External cursor/foreground change before native hover move");
  if(!SetCursorPos(target.X,target.Y)) throw new Exception("Owned hover cursor move failed");
  POINT after;if(!GetCursorPos(out after)||after.X!=target.X||after.Y!=target.Y) throw new Exception("Owned hover cursor did not reach target");
  OwnedPoint(h,pid,after);return new int[]{after.X,after.Y};
 }
 public static string Restore(int ox,int oy,int lx,int ly,int tx,int ty) {
  POINT p;if(!GetCursorPos(out p)) throw new Exception("Cursor unavailable for restoration");
  if(p.X==ox&&p.Y==oy) return "already-original";
  if(!((p.X==lx&&p.Y==ly)||(p.X==tx&&p.Y==ty))) return "external-position-preserved";
  if(!SetCursorPos(ox,oy)) throw new Exception("Cursor restoration failed");
  if(!GetCursorPos(out p)||p.X!=ox||p.Y!=oy) throw new Exception("Cursor restoration did not remain at original point");
  return "restored";
 }
}
