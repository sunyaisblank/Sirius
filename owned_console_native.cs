using System;
using System.Text;
using System.Runtime.InteropServices;
using System.ComponentModel;

// Task-local diagnostic adapter for an unchanged console executable.
// Screen cells are preserved separately; they are not original stdout bytes.
public static class SiriusOwnedConsole {
 [StructLayout(LayoutKind.Sequential, CharSet=CharSet.Unicode)]
 struct StartupInfo {
  public uint cb; public string reserved, desktop, title;
  public uint x,y,xSize,ySize,xCountChars,yCountChars,fillAttribute,flags;
  public ushort showWindow,reserved2; public IntPtr reservedBytes,input,output,error;
 }
 [StructLayout(LayoutKind.Sequential)]
 struct ProcessInfo { public IntPtr process,thread; public uint pid,tid; }
 [StructLayout(LayoutKind.Sequential)]
 public struct Coord { public short x,y; public Coord(short x,short y){this.x=x;this.y=y;} }
 [StructLayout(LayoutKind.Sequential)]
 struct SmallRect { public short left,top,right,bottom; }
 [StructLayout(LayoutKind.Sequential)]
 struct BufferInfo { public Coord size,cursor; public ushort attributes; public SmallRect window; public Coord maximum; }
 [DllImport("kernel32.dll",CharSet=CharSet.Unicode,SetLastError=true)]
 static extern bool CreateProcessW(string app,StringBuilder command,IntPtr pa,IntPtr ta,bool inherit,uint flags,IntPtr env,string cwd,ref StartupInfo si,out ProcessInfo pi);
 [DllImport("kernel32.dll",SetLastError=true)] static extern uint ResumeThread(IntPtr thread);
 [DllImport("kernel32.dll",SetLastError=true)] static extern bool CloseHandle(IntPtr h);
 [DllImport("kernel32.dll",SetLastError=true)] static extern bool FreeConsole();
 [DllImport("kernel32.dll",SetLastError=true)] static extern bool AttachConsole(uint pid);
 [DllImport("kernel32.dll",SetLastError=true)] static extern uint GetConsoleProcessList([Out]uint[] ids,uint count);
 [DllImport("kernel32.dll",CharSet=CharSet.Unicode,SetLastError=true)]
 static extern IntPtr CreateFileW(string path,uint access,uint share,IntPtr sa,uint creation,uint flags,IntPtr template);
 [DllImport("kernel32.dll",SetLastError=true)] static extern bool GetConsoleScreenBufferInfo(IntPtr h,out BufferInfo info);
 [DllImport("kernel32.dll",CharSet=CharSet.Unicode,SetLastError=true)]
 static extern bool ReadConsoleOutputCharacterW(IntPtr h,[Out]StringBuilder text,uint length,Coord origin,out uint read);
 static ProcessInfo launch;
 static IntPtr screen=IntPtr.Zero;
 static bool attached=false;
 public static int LastAttachError {get;private set;}
 static int minimumWidth,observedWidth,observedHeight;
 public static void ConfigureWidth(int width) {
  if(width<1 || width>4096 || launch.process!=IntPtr.Zero || attached)throw new Exception("Invalid prelaunch width contract");
  minimumWidth=width;
 }
 static Exception Error(){return new Win32Exception(Marshal.GetLastWin32Error());}
 public static uint StartSuspended(string app,string arguments,string cwd) {
  if(launch.process!=IntPtr.Zero)throw new Exception("Console child already created");
  var si=new StartupInfo();si.cb=(uint)Marshal.SizeOf(typeof(StartupInfo));si.flags=1;si.showWindow=0;
  // Real console, suspended until the exact process birth has been recorded.
  if(!CreateProcessW(app,new StringBuilder("\""+app+"\" "+arguments),IntPtr.Zero,IntPtr.Zero,false,0x14,IntPtr.Zero,cwd,ref si,out launch))throw Error();
  return launch.pid;
 }
 public static void Resume() {
  if(launch.thread==IntPtr.Zero || ResumeThread(launch.thread)!=1)throw Error();
  if(!CloseHandle(launch.thread))throw Error();launch.thread=IntPtr.Zero;
 }
 public static bool TryAttach(uint pid) {
  if(attached)throw new Exception("Console observer already attached");
  if(!FreeConsole())throw Error();
  if(!AttachConsole(pid)){LastAttachError=Marshal.GetLastWin32Error();return false;}
  LastAttachError=0;
  attached=true;
  screen=CreateFileW("CONOUT$",0x80000000,3,IntPtr.Zero,3,0,IntPtr.Zero);
  if(screen==new IntPtr(-1)) {screen=IntPtr.Zero;throw Error();}
  // Read-only screen access; actual geometry is bound on the first snapshot.
  return true;
 }
 public sealed class Snapshot {
  public int width,height,cursor_x,cursor_y,rows_read;
  public uint[] attached_pids;
  public string cells,text;
 }
 public static Snapshot Read(uint owner,uint child) {
  if(!attached || screen==IntPtr.Zero)throw new Exception("No owned console attachment");
  var ids=new uint[16];uint count=GetConsoleProcessList(ids,(uint)ids.Length);
  if(count==0 || count>ids.Length)throw new Exception("Console membership unavailable or oversized");
  bool hasOwner=false;for(int i=0;i<count;i++){if(ids[i]==owner)hasOwner=true;else if(ids[i]!=child)throw new Exception("Unowned process shares console");}
  if(!hasOwner)throw new Exception("Console observer membership absent");
  BufferInfo info;if(!GetConsoleScreenBufferInfo(screen,out info))throw Error();
  if(minimumWidth<1 || info.size.x<minimumWidth || info.size.x>4096 || info.size.y<2 || info.cursor.y<0 || info.cursor.y>=info.size.y-1 || info.size.x*(info.cursor.y+1)>2097152)throw new Exception("Console geometry/read/wrap/scroll bound invalid: "+info.size.x+"x"+info.size.y);
  if(observedWidth==0){observedWidth=info.size.x;observedHeight=info.size.y;}
  if(info.size.x!=observedWidth || info.size.y!=observedHeight)throw new Exception("Console geometry changed during owned session");
  uint length=(uint)(info.size.x*(info.cursor.y+1));var cells=new StringBuilder((int)length);uint read;
  if(!ReadConsoleOutputCharacterW(screen,cells,length,new Coord(0,0),out read) || read!=length)throw Error();
  string raw=cells.ToString();var text=new StringBuilder();
  for(int row=0;row<=info.cursor.y;row++) {
   string line=raw.Substring(row*info.size.x,info.size.x);
   if(line[line.Length-1]!=' ')throw new Exception("Console line wrapping ambiguity");
   text.AppendLine(line.TrimEnd(' '));
  }
  Array.Resize(ref ids,(int)count);
  return new Snapshot{width=info.size.x,height=info.size.y,cursor_x=info.cursor.x,cursor_y=info.cursor.y,rows_read=info.cursor.y+1,attached_pids=ids,cells=raw,text=text.ToString()};
 }
 public static void Detach() {
  Exception error=null;
  if(screen!=IntPtr.Zero){if(!CloseHandle(screen))error=Error();screen=IntPtr.Zero;}
  if(attached){if(!FreeConsole() && error==null)error=Error();attached=false;}
  if(launch.thread!=IntPtr.Zero){if(!CloseHandle(launch.thread) && error==null)error=Error();launch.thread=IntPtr.Zero;}
  if(launch.process!=IntPtr.Zero){if(!CloseHandle(launch.process) && error==null)error=Error();launch.process=IntPtr.Zero;}
  if(error!=null)throw error;
 }
}
