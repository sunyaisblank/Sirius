param([Parameter(Mandatory=$true)][string]$Tools,[Parameter(Mandatory=$true)][string]$Output,[Parameter(Mandatory=$true)][ValidateSet('profiler-trace-owned-click')][string]$Action,[Parameter(Mandatory=$true)][string]$ExpectedHash,[Parameter(Mandatory=$true)][string]$LiveRevision,[Parameter(Mandatory=$true)][string]$Trace,[Parameter(Mandatory=$true)][string]$TraceHash)
$ErrorActionPreference='Stop'
if (Test-Path -LiteralPath $Output) { throw 'One-shot tool output already exists' }
$Output=(New-Item -ItemType Directory -Path $Output).FullName
$tempFull=(New-Item -ItemType Directory -Path (Join-Path $Output 'temp')).FullName
$env:TEMP=$tempFull;$env:TMP=$tempFull
$tool='RadeonGPUProfiler.exe'
$exe=Join-Path $Tools $tool
$ownerProcess=Get-Process -Id $PID;$ownerHandle=$ownerProcess.Handle
$ownerStart=$ownerProcess.StartTime.ToUniversalTime().ToString('o')
$bootstrap=[ordered]@{owner_pid=$PID;owner_start_time_utc=$ownerStart;expected_test_executable=$exe;scope='Early exact GUI observer birth; no service/capture/scientific launch authorized.'}
$bootstrap | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $tempFull 'bootstrap.json') -Encoding UTF8
Add-Type -AssemblyName UIAutomationClient,UIAutomationTypes
if (!(Test-Path -LiteralPath $Trace) -or (Get-Item -LiteralPath $Trace).Length -le 0 -or (Get-Item -LiteralPath $Trace).Length -gt 268435456 -or (Get-FileHash -LiteralPath $Trace -Algorithm SHA256).Hash.ToLowerInvariant() -ne $TraceHash) {throw 'Pinned trace input mismatch'}
Add-Type -ReferencedAssemblies System.Drawing,System,System.Core -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Text;
public static class SiriusOwnedExit {
 [DllImport("kernel32.dll", SetLastError=true)]
 [return: MarshalAs(UnmanagedType.Bool)]
 private static extern bool GetExitCodeProcess(IntPtr hProcess, out uint code);
 public static uint Read(IntPtr handle) {
  uint code;
  if (!GetExitCodeProcess(handle, out code)) throw new System.ComponentModel.Win32Exception(Marshal.GetLastWin32Error());
  if (code == 259) throw new Exception("Owned process is still active");
  return code;
 }
}
public sealed class SiriusWindowInfo {
 public long handle; public uint pid; public string title; public int width; public int height;
 public bool print_window_succeeded; public string image_path; public long image_bytes; public string image_sha256;
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
 [StructLayout(LayoutKind.Sequential)] private struct POINT {public int X,Y;}
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
 public static SiriusWindowInfo[] Capture(uint pid, string output) {
  var rows = new List<SiriusWindowInfo>(); Exception failure=null;
  bool enumerated = EnumWindows(delegate(IntPtr h, IntPtr data) {
   try {
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
     try {
      GetWindowThreadProcessId(h,out actual);
      if(actual!=pid) throw new Exception("Window owner changed before PrintWindow");
      row.print_window_succeeded=PrintWindow(h,dc,2);
      GetWindowThreadProcessId(h,out actual);
      if(actual!=pid) throw new Exception("Window owner changed during PrintWindow");
     } finally {graphics.ReleaseHdc(dc);}
    }
    if(row.print_window_succeeded) {
     row.image_path=Path.Combine(output,"owned-window-"+h.ToInt64().ToString("x")+".png");
     bitmap.Save(row.image_path,ImageFormat.Png);
     row.image_bytes=new FileInfo(row.image_path).Length;
     using(var algorithm=System.Security.Cryptography.SHA256.Create())
      using(var stream=File.OpenRead(row.image_path))
       row.image_sha256=BitConverter.ToString(algorithm.ComputeHash(stream)).Replace("-","").ToLowerInvariant();
    }
   }
   rows.Add(row); return true;
   } catch(Exception e) {failure=e;return false;}
  },IntPtr.Zero);
  if(failure!=null) throw failure;
  if(!enumerated) throw new Exception("Owned window enumeration failed");
  return rows.ToArray();
 }
}

'@
$existing=@(Get-Process | Where-Object {$_.ProcessName -match '^Radeon(Developer|GPUProfiler)'})
if ($existing.Count) { throw 'Existing developer/profiler process: tool help declined' }
$hash=(Get-FileHash -LiteralPath $exe -Algorithm SHA256).Hash.ToLowerInvariant()
if ($hash -ne $ExpectedHash -or $LiveRevision -notmatch '^[0-9a-f]{40}$') { throw 'Pinned tool/source mismatch' }
$signature=Get-AuthenticodeSignature -LiteralPath $exe
$sig=[ordered]@{status=$signature.Status.ToString();subject=$signature.SignerCertificate.Subject;thumbprint=$signature.SignerCertificate.Thumbprint;file_version=(Get-Item -LiteralPath $exe).VersionInfo.FileVersion;product_version=(Get-Item -LiteralPath $exe).VersionInfo.ProductVersion}
$sig | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'signature.json') -Encoding UTF8
if ($signature.Status -ne 'Valid') { throw 'Official tool signature did not verify' }
$record=[ordered]@{action=$Action;owner_pid=$PID;owner_start_time_utc=$ownerStart;test_executable=$exe;test_executable_sha256=$hash;arguments=@($Trace);trace_sha256=$TraceHash;trace_path=$Trace;intentional_owned_stop=$false;window_observation=@();live_source_revision=$LiveRevision;polling_guard_seconds=60;blocked_call_bridge_guard_seconds=120;rss_guard_bytes=1073741824;outer_stop=$false;status='running';scope='Owned saved-profile tabs and one bounded HWND client mouse-message request/result observation only; no global input API or assumed effect; positional file opening is a bounded probe of the documented executable-drop route. HWND PID checked before and after PrintWindow; synchronous call is subject to outer bridge recovery. Deliberate owned stop after sealed PNG observation; no service/capture/scientific result or clock/driver changes.'}
$process=$null;$watch=[Diagnostics.Stopwatch]::StartNew();$peak=0;$errors=New-Object 'System.Collections.Generic.List[string]'
try {
 $process=Start-Process -FilePath $exe -ArgumentList @(('"'+$Trace+'"')) -WorkingDirectory $Tools -NoNewWindow -RedirectStandardOutput (Join-Path $Output 'stdout.log') -RedirectStandardError (Join-Path $Output 'stderr.log') -PassThru
 $nativeHandle=$process.Handle
 $record.ui_unavailable_pid_zero=0
 $record.test_pid=$process.Id;$record.test_start_time_utc=$process.StartTime.ToUniversalTime().ToString('o')
 $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 while (!$process.HasExited) {
  $process.Refresh();if ($process.HasExited) { break }
  $peak=[Math]::Max($peak,$process.WorkingSet64)
  if ($watch.Elapsed.TotalSeconds -gt 60 -or $process.WorkingSet64 -gt 1073741824) {
   $record.outer_stop=$true;$process.Kill()
   if (!$process.WaitForExit(5000)) { throw 'Owned tool did not exit after guard stop' };break
  }
  if ($watch.Elapsed.TotalSeconds -ge 5) {
   $rows=@([SiriusWindow]::Capture([uint32]$process.Id,$Output))
   if ($rows.Count) {
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
    if ($rows.Count -ne 1 -or $process.HasExited) {throw 'Navigation requires exactly one live owned root'}
    $root=[System.Windows.Automation.AutomationElement]::FromHandle([IntPtr]$rows[0].handle)
    if ($root.Current.ProcessId -ne $process.Id) {throw 'Navigation root owner mismatch'}
    $record.navigation_history=@()
    foreach ($tabLabel in @('EVENTS','Event timing')) {
    $tabName=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::NameProperty,$tabLabel)
    $tabType=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::ControlTypeProperty,[System.Windows.Automation.ControlType]::TabItem)
    $tabCondition=New-Object System.Windows.Automation.AndCondition($tabName,$tabType)
    $tabs=$root.FindAll([System.Windows.Automation.TreeScope]::Descendants,$tabCondition)
    if ($tabs.Count -ne 1) {throw 'Expected unique observed EVENTS tab'}
    $tab=$tabs[0];$c=$tab.Current
    if ($process.HasExited -or $c.ProcessId -ne $process.Id -or !$c.IsEnabled -or $c.IsOffscreen) {throw 'EVENTS tab owner/lifetime/state mismatch'}
    $selection=$null
    if (!$tab.TryGetCurrentPattern([System.Windows.Automation.SelectionItemPattern]::Pattern,[ref]$selection)) {throw 'EVENTS selection pattern unavailable'}
    $record.navigation_requested=$tabLabel+' SelectionItemPattern.Select on exact owned root'
    $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
    if ($process.HasExited -or $root.Current.ProcessId -ne $process.Id -or $tab.Current.ProcessId -ne $process.Id) {throw 'Navigation ownership changed before select'}
    $selection.Select()
    Start-Sleep -Milliseconds 500
    if ($process.HasExited -or $root.Current.ProcessId -ne $process.Id -or $tab.Current.ProcessId -ne $process.Id -or !$selection.Current.IsSelected) {throw 'EVENTS selection did not remain owned/selected'}
    $record.navigation_history+=@($tabLabel)
    }
    $record.navigation_completed=$true


    $record.event_observations=@()
    foreach ($eventName in @('4 vkCmdDispatch(3, 1, 1)')) {
     $label=if ($eventName -eq '4 vkCmdDispatch(3, 1, 1)') {'transport'} else {'endpoint'}
     if ($process.HasExited -or $root.Current.ProcessId -ne $process.Id) {throw 'Event selection root owner/lifetime mismatch'}
     $eventNameCondition=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::NameProperty,$eventName)
     $eventTypeCondition=New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::ControlTypeProperty,[System.Windows.Automation.ControlType]::TreeItem)
     $eventCondition=New-Object System.Windows.Automation.AndCondition($eventNameCondition,$eventTypeCondition)
     $events=$root.FindAll([System.Windows.Automation.TreeScope]::Descendants,$eventCondition)
     if ($events.Count -ne 1) {throw 'Observed dispatch must be unique'}
     $event=$events[0];$c=$event.Current
     if ($process.HasExited -or $c.ProcessId -ne $process.Id -or !$c.IsEnabled -or $c.IsOffscreen) {throw 'Dispatch selection owner/lifetime/state mismatch'}
     $record.event_selection_requested=$eventName
     $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
     if ($process.HasExited -or $event.Current.ProcessId -ne $process.Id -or $root.Current.ProcessId -ne $process.Id) {throw 'Dispatch selection ownership changed'}
     $r=$event.Current.BoundingRectangle
     $record.owned_message_target=[ordered]@{hwnd=$rows[0].handle;pid=$process.Id;event_name=$eventName;rect=@($r.X,$r.Y,$r.Width,$r.Height)}
     $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
     if ($process.HasExited -or $event.Current.ProcessId -ne $process.Id -or $root.Current.ProcessId -ne $process.Id) {throw 'Owned mouse-message target lifetime changed'}
     $record.owned_mouse_message_result=[SiriusWindow]::ClickOwned($rows[0].handle,[uint32]$process.Id,$r.X,$r.Y,$r.Width,$r.Height)
     Start-Sleep -Milliseconds 500
     if ($process.HasExited -or $root.Current.ProcessId -ne $process.Id) {throw 'Owned root/lifetime changed after selection request'}
     $record.selection_result=[ordered]@{original_event_pid=$event.Current.ProcessId;requested_action='Owned HWND client WM_LBUTTONDOWN/UP'}
     $fresh=$root.FindAll([System.Windows.Automation.TreeScope]::Descendants,$eventCondition)
     $record.selection_result.fresh_count=$fresh.Count
     if ($fresh.Count -ne 1 -or $process.HasExited -or $fresh[0].Current.ProcessId -ne $process.Id) {throw 'Fresh dispatch reference is not uniquely owned'}
     $freshSelection=$null
     if (!$fresh[0].TryGetCurrentPattern([System.Windows.Automation.SelectionItemPattern]::Pattern,[ref]$freshSelection)) {throw 'Fresh dispatch selection pattern unavailable'}
     $record.selection_result.fresh_is_selected=$freshSelection.Current.IsSelected
     $record.selection_result.available_patterns=@($fresh[0].GetSupportedPatterns() | ForEach-Object {$_.ProgrammaticName})
     if ($process.HasExited -or $root.Current.ProcessId -ne $process.Id -or $fresh[0].Current.ProcessId -ne $process.Id) {throw 'Fresh selection observation ownership changed'}

    if ($process.HasExited) {throw 'Owned profiler exited during window observation'}
    $record.current_window_observation=$rows
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
      if ($process.HasExited) {throw 'Owned profiler exited during ready UI observation'}
      if ($c.ProcessId -eq 0) {$record.ui_unavailable_pid_zero++;continue}
      if ($c.ProcessId -ne $process.Id) {throw 'UI child foreign owner mismatch'}
      $value=$null;$pattern=$null
      if ($element.TryGetCurrentPattern([System.Windows.Automation.ValuePattern]::Pattern,[ref]$pattern)) {$value=$pattern.Current.Value}
      if ($c.Name.Length -gt 4096 -or ($null -ne $value -and $value.Length -gt 4096)) {throw 'Owned UI text exceeds diagnostic bound'}
      $automation.Add([pscustomobject]@{name=$c.Name;value=$value;control=$c.ControlType.ProgrammaticName;automation_id=$c.AutomationId;class=$c.ClassName;enabled=$c.IsEnabled;offscreen=$c.IsOffscreen;rect=@($c.BoundingRectangle.X,$c.BoundingRectangle.Y,$c.BoundingRectangle.Width,$c.BoundingRectangle.Height)})
     }
    }
    $automation.ToArray() | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output ($label+'-automation.json')) -Encoding UTF8
    $record.owned_ui_elements=$automation.Count
    $main=@($automation.ToArray() | Where-Object {$_.automation_id -eq 'MainWindow.main_widget_.main_tab_widget_'})
    $loading=@($automation.ToArray() | Where-Object {$_.automation_id -eq 'MainWindow.file_loading_widget_' -and !$_.offscreen})
    $record.loaded_profile_observed=($main.Count -eq 1 -and $main[0].enabled -and $loading.Count -eq 0)
    if (!$record.loaded_profile_observed) {throw 'Loaded profile readiness lost after Event timing selection'}
    # The first screenshot may precede the loader transition; recapture now.
    $rows=@([SiriusWindow]::Capture([uint32]$process.Id,$Output))
    if ($process.HasExited) {throw 'Owned profiler exited during loaded-profile observation'}
    $record.current_window_observation=$rows
    if (!@($rows | Where-Object {$_.print_window_succeeded -and $_.image_bytes -gt 0 -and $_.image_sha256}).Count) {throw 'No successfully saved owned-window PNG'}
    foreach ($row in $rows) {
     if ($row.print_window_succeeded -and $row.image_path) {
      $renamed=Join-Path $Output ($label+'-'+[IO.Path]::GetFileName($row.image_path))
      if (Test-Path -LiteralPath $renamed) {throw 'One-shot event image already exists'}
      Move-Item -LiteralPath $row.image_path -Destination $renamed
      $row.image_path=$renamed
     }
    }
    $record.event_observations+=@([pscustomobject]@{event_name=$eventName;label=$label;owned_ui_elements=$automation.Count;automation_path=(Join-Path $Output ($label+'-automation.json'));window_observation=$rows})
    }
    $record.window_observation=@($record.event_observations | ForEach-Object {$_.window_observation})
    $record.intentional_owned_stop=$true
    $process.Kill();if (!$process.WaitForExit(5000)) {throw 'Observed owned profiler did not terminate'};break
   }
  }
  Start-Sleep -Milliseconds 50
 }
 if (!$process.WaitForExit(5000)) { throw 'Owned tool did not complete its exit' }
 $record.returncode=[SiriusOwnedExit]::Read($nativeHandle)
 $record.exit_code_source='GetExitCodeProcess on preacquired owned handle after WaitForExit'
} catch { $record.exception=$_.Exception.Message;$record.status='failed' }
finally {
 $record.owned_process_absent=($null -eq $process)
 $record.tool_already_exited_on_finally_entry=$null;$record.exception_cleanup_kill_requested=$false
 if ($process) {
  try {$record.tool_already_exited_on_finally_entry=$process.HasExited;if (!$record.tool_already_exited_on_finally_entry) {$record.exception_cleanup_kill_requested=$true;$process.Kill();if (!$process.WaitForExit(5000)) {throw 'Owned help tool did not terminate'}} } catch {$errors.Add($_.Exception.Message)}
  try {$record.owned_process_absent=$process.HasExited;if ($record.owned_process_absent) {$record.terminal_returncode=[SiriusOwnedExit]::Read($nativeHandle);$record.terminal_exit_code_source='Held native handle after terminal wait, including failed observation'}} catch {$errors.Add($_.Exception.Message)}
 }
 $record.elapsed_seconds=$watch.Elapsed.TotalSeconds;$record.peak_sampled_rss_bytes=$peak
 $record.trace_unchanged=$false
 try {$record.trace_unchanged=((Get-FileHash -LiteralPath $Trace -Algorithm SHA256).Hash.ToLowerInvariant() -eq $TraceHash)} catch {$errors.Add('trace-posthash: '+$_.Exception.Message)}
 $record.tool_unchanged=$false
 try {$record.tool_unchanged=((Get-FileHash -LiteralPath $exe -Algorithm SHA256).Hash.ToLowerInvariant() -eq $hash)} catch {$errors.Add('posthash: '+$_.Exception.Message)}
 $record.accessible_remaining_developer_processes=@();$record.developer_census_completed=$false
 try {
  $record.accessible_remaining_developer_processes=@(Get-Process | Where-Object {$_.ProcessName -match '^Radeon(Developer|GPUProfiler)'} | ForEach-Object {@{pid=$_.Id;name=$_.ProcessName;start_utc=$_.StartTime.ToUniversalTime().ToString('o')}})
  $record.developer_census_completed=$true
 } catch {$errors.Add('developer-census: '+$_.Exception.Message)}
 $record.cleanup_errors=$errors.ToArray()
 if ($record.status -ne 'failed') {$record.status=if ($record.outer_stop) {'invalidated'} elseif ($record.intentional_owned_stop -and $record.window_observation.Count) {'observed-and-stopped'} else {'no-window-observation'}}
 if (!$record.trace_unchanged -or $errors.Count -or !$record.owned_process_absent -or !$record.tool_unchanged -or !$record.developer_census_completed -or $record.accessible_remaining_developer_processes.Count) {$record.status='failed'}
 try {$record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8} catch {
  $record.status='failed';$record.owner_receipt_error=$_.Exception.Message
  [Console]::Error.WriteLine('Terminal tool owner receipt failed: '+$record.owner_receipt_error)
 }
}
$record | ConvertTo-Json -Depth 8
if ($record.event_observations.Count -ne 1 -or !$record.navigation_completed -or !$record.loaded_profile_observed -or $record.status -ne 'observed-and-stopped' -or !$record.intentional_owned_stop -or !$record.window_observation.Count) {exit 1}
