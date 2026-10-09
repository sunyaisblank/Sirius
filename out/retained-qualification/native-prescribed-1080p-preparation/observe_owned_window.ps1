param(
    [Parameter(Mandatory=$true)][int]$ViewerProcessId,
    [Parameter(Mandatory=$true)][long]$StartTimeTicks,
    [Parameter(Mandatory=$true)][string]$ExecutableSha256,
    [Parameter(Mandatory=$true)][ValidateSet('window','capture','inputs','escape')][string]$Action,
    [Parameter(Mandatory=$true)][string]$RecordPath,
    [string]$CapturePath
)
# Diagnostic preparation only. The existing collector creates/owns the viewer;
# this helper neither starts nor kills it and never uses the global input stream.
$ErrorActionPreference = 'Stop'
[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new($false)
$utf8 = [System.Text.UTF8Encoding]::new($false)
if (Test-Path -LiteralPath $RecordPath) { throw 'Refuse existing observation record' }
if ($Action -eq 'capture' -and ([string]::IsNullOrEmpty($CapturePath) -or (Test-Path -LiteralPath $CapturePath))) { throw 'Capture needs a new path' }
Add-Type -TypeDefinition @'
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
public static class SiriusOwnedWindow {
    public delegate bool EnumProc(IntPtr window, IntPtr parameter);
    [StructLayout(LayoutKind.Sequential)] public struct Point { public int x, y; }
    [StructLayout(LayoutKind.Sequential)] public struct Rect { public int left, top, right, bottom; }
    [DllImport("user32.dll")] public static extern bool EnumWindows(EnumProc callback, IntPtr parameter);
    [DllImport("user32.dll")] public static extern uint GetWindowThreadProcessId(IntPtr window, out uint process);
    [DllImport("user32.dll")] public static extern bool IsWindowVisible(IntPtr window);
    [DllImport("user32.dll")] public static extern bool IsIconic(IntPtr window);
    [DllImport("user32.dll")] public static extern bool GetClientRect(IntPtr window, out Rect rectangle);
    [DllImport("user32.dll")] public static extern bool ClientToScreen(IntPtr window, ref Point point);
    [DllImport("user32.dll")] public static extern IntPtr GetForegroundWindow();
    [DllImport("user32.dll")] public static extern IntPtr WindowFromPoint(Point point);
    [DllImport("user32.dll")] public static extern bool SetProcessDPIAware();
    [DllImport("user32.dll", SetLastError=true)] public static extern bool PostMessageW(IntPtr window, uint message, UIntPtr wparam, IntPtr lparam);
    public static IntPtr[] VisibleWindows(uint process) {
        var found = new List<IntPtr>();
        if (!EnumWindows(delegate(IntPtr w, IntPtr p) {
            uint owner; GetWindowThreadProcessId(w, out owner);
            if (owner == process && IsWindowVisible(w)) found.Add(w);
            return true;
        }, IntPtr.Zero)) throw new InvalidOperationException("EnumWindows failed");
        return found.ToArray();
    }
}
'@
$p = Get-Process -Id $ViewerProcessId
$null = $p.Handle
function Assert-Owner {
    $p.Refresh()
    if ($p.HasExited -or $p.StartTime.ToUniversalTime().Ticks -ne $StartTimeTicks) { throw 'Owned process lifetime changed' }
    if ((Get-FileHash -Algorithm SHA256 -LiteralPath $p.MainModule.FileName).Hash.ToLowerInvariant() -ne $ExecutableSha256.ToLowerInvariant()) { throw 'Owned executable changed' }
}
Assert-Owner
$windows = @([SiriusOwnedWindow]::VisibleWindows([uint32]$ViewerProcessId))
if ($windows.Count -ne 1) { throw 'Expected exactly one visible top-level window belonging to the owned viewer' }
$window = $windows[0]
function Assert-Window {
    Assert-Owner
    [uint32]$owner = 0
    $null = [SiriusOwnedWindow]::GetWindowThreadProcessId($window, [ref]$owner)
    if ($owner -ne $ViewerProcessId -or -not [SiriusOwnedWindow]::IsWindowVisible($window)) { throw 'Owned HWND changed or became invisible' }
}
function Post-Owned([uint32]$message, [uint64]$word, [long]$bits) {
    Assert-Window
    if (-not [SiriusOwnedWindow]::PostMessageW($window, $message, [UIntPtr]::new($word), [IntPtr]::new($bits))) { throw ('PostMessage failed: ' + [Runtime.InteropServices.Marshal]::GetLastWin32Error()) }
}
$stateLines = @(& quser.exe $p.SessionId 2>&1 | ForEach-Object { [string]$_ })
$sessionState = 'unparsed'
if ($stateLines.Count -gt 1 -and $stateLines[1] -match '\s(Active|Disc|Disconnected|Connect|Conn|Listen)\s') { $sessionState = $Matches[1] }
$record = [ordered]@{schema_version=1; pid=$ViewerProcessId; start_time_ticks=$StartTimeTicks; executable=$p.MainModule.FileName; executable_sha256=$ExecutableSha256.ToLowerInvariant(); hwnd=$window.ToInt64(); session_id=$p.SessionId; session_state=$sessionState; action=$Action; started_utc=[DateTime]::UtcNow.ToString('o'); input_scope='owned HWND PostMessage; no SendInput/hardware-input claim'}
if ($Action -eq 'capture') {
    if ($sessionState -ne 'Active' -or [SiriusOwnedWindow]::IsIconic($window)) { throw 'Visible capture requires an active desktop and a non-minimized viewer' }
    $null = [SiriusOwnedWindow]::SetProcessDPIAware()
    Assert-Window
    if ([SiriusOwnedWindow]::GetForegroundWindow() -ne $window) { throw 'Capture requires owned viewer already foreground; helper does not steal focus' }
    $rect = [SiriusOwnedWindow+Rect]::new()
    $origin = [SiriusOwnedWindow+Point]::new()
    if (-not [SiriusOwnedWindow]::GetClientRect($window, [ref]$rect) -or -not [SiriusOwnedWindow]::ClientToScreen($window, [ref]$origin)) { throw 'Client bounds unavailable' }
    $width = $rect.right - $rect.left; $height = $rect.bottom - $rect.top
    if ($width -le 0 -or $height -le 0 -or $width -gt 8192 -or $height -gt 8192) { throw 'Invalid finite capture bounds' }
    # These samples reject ordinary occlusion; they are not a proof that every
    # pixel stayed unobscured between the checks and desktop capture.
    foreach ($xy in @(@(1,1),@(($width-2),1),@(1,($height-2)),@(($width-2),($height-2)),@([int]($width/2),[int]($height/2)))) {
        $point = [SiriusOwnedWindow+Point]::new(); $point.x=$origin.x+$xy[0]; $point.y=$origin.y+$xy[1]
        if ([SiriusOwnedWindow]::WindowFromPoint($point) -ne $window) { throw 'Owned client is occluded at a sampled point' }
    }
    Add-Type -AssemblyName System.Drawing
    $bitmap = [Drawing.Bitmap]::new($width,$height)
    $graphics = [Drawing.Graphics]::FromImage($bitmap)
    try {
        $graphics.CopyFromScreen($origin.x,$origin.y,0,0,[Drawing.Size]::new($width,$height))
        Assert-Window
        $bitmap.Save($CapturePath,[Drawing.Imaging.ImageFormat]::Png)
    } finally { $graphics.Dispose(); $bitmap.Dispose() }
    $record['client_screen_rectangle']=@($origin.x,$origin.y,$width,$height)
    $record['capture_path']=$CapturePath
    $record['capture_sha256']=(Get-FileHash -Algorithm SHA256 -LiteralPath $CapturePath).Hash.ToLowerInvariant()
    $record['capture_scope']='desktop-composited owned client pixels; not a GL upload/swap timestamp or scientific pixel reference'
} elseif ($Action -eq 'inputs') {
    if ($sessionState -ne 'Active') { throw 'Window-input observation requires active viewer desktop' }
    # GLFW maps virtual keys when scan code is zero. Bit31 determines release.
    Post-Owned 0x0100 0x57 1
    try { Start-Sleep -Milliseconds 200 } finally { Post-Owned 0x0101 0x57 3221225473 }
    Post-Owned 0x0200 0 ((120 -shl 16) -bor 120)
    Post-Owned 0x0201 1 ((120 -shl 16) -bor 120)
    Post-Owned 0x0200 1 ((160 -shl 16) -bor 180)
    Post-Owned 0x0202 0 ((160 -shl 16) -bor 180)
    Post-Owned 0x020A (120 -shl 16) 0
    $record['posted_events']=@('W down/up','cursor120,120','left down','drag180,160','left up','wheel+120')
} elseif ($Action -eq 'escape') {
    Post-Owned 0x0100 0x1B 1
    # Press alone is the production exit trigger; posting a release after it
    # could race a legitimately completed/closed owned process.
    $record['posted_events']=@('Escape down')
    $record['termination_scope']='request only; existing owner must observe normal exit and WaitForExit, otherwise owned kill fallback is not a cooperative-cancel pass'
}
if ($Action -ne 'escape') { Assert-Owner }
$record['completed_utc']=[DateTime]::UtcNow.ToString('o')
[IO.File]::WriteAllText($RecordPath,($record | ConvertTo-Json -Depth 6),$utf8)
$p.Dispose()
