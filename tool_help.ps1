param([Parameter(Mandatory=$true)][string]$Tools,[Parameter(Mandatory=$true)][string]$Output,[Parameter(Mandatory=$true)][ValidateSet('panel-help','panel-version','profiler-help')][string]$Action,[Parameter(Mandatory=$true)][string]$ExpectedHash,[Parameter(Mandatory=$true)][string]$LiveRevision)
$ErrorActionPreference='Stop'
if (Test-Path -LiteralPath $Output) { throw 'One-shot tool output already exists' }
$Output=(New-Item -ItemType Directory -Path $Output).FullName
$tempFull=(New-Item -ItemType Directory -Path (Join-Path $Output 'temp')).FullName
$env:TEMP=$tempFull;$env:TMP=$tempFull
$tool=if ($Action -eq 'profiler-help') {'RadeonGPUProfiler.exe'} else {'RadeonDeveloperPanelCLI.exe'}
$argument=if ($Action -eq 'panel-version') {'--version'} else {'--help'}
$exe=Join-Path $Tools $tool
$ownerProcess=Get-Process -Id $PID;$ownerHandle=$ownerProcess.Handle
$ownerStart=$ownerProcess.StartTime.ToUniversalTime().ToString('o')
$bootstrap=[ordered]@{owner_pid=$PID;owner_start_time_utc=$ownerStart;expected_test_executable=$exe;scope='Early exact tool-help observer birth; no service/capture/scientific launch authorized.'}
$bootstrap | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $tempFull 'bootstrap.json') -Encoding UTF8
Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
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
'@
$existing=@(Get-Process | Where-Object {$_.ProcessName -match '^Radeon(Developer|GPUProfiler)'})
if ($existing.Count) { throw 'Existing developer/profiler process: tool help declined' }
$hash=(Get-FileHash -LiteralPath $exe -Algorithm SHA256).Hash.ToLowerInvariant()
if ($hash -ne $ExpectedHash -or $LiveRevision -notmatch '^[0-9a-f]{40}$') { throw 'Pinned tool/source mismatch' }
$signature=Get-AuthenticodeSignature -LiteralPath $exe
$sig=[ordered]@{status=$signature.Status.ToString();subject=$signature.SignerCertificate.Subject;thumbprint=$signature.SignerCertificate.Thumbprint;file_version=(Get-Item -LiteralPath $exe).VersionInfo.FileVersion;product_version=(Get-Item -LiteralPath $exe).VersionInfo.ProductVersion}
$sig | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'signature.json') -Encoding UTF8
if ($signature.Status -ne 'Valid') { throw 'Official tool signature did not verify' }
$record=[ordered]@{action=$Action;owner_pid=$PID;owner_start_time_utc=$ownerStart;test_executable=$exe;test_executable_sha256=$hash;arguments=@($argument);live_source_revision=$LiveRevision;outer_guard_seconds=30;rss_guard_bytes=1073741824;outer_stop=$false;status='running';scope='Pinned official executable help/version only; no service/capture/scientific result or developer-mode/clock/driver experiment authorized.'}
$process=$null;$watch=[Diagnostics.Stopwatch]::StartNew();$peak=0;$errors=New-Object 'System.Collections.Generic.List[string]'
try {
 $process=Start-Process -FilePath $exe -ArgumentList $argument -WorkingDirectory $Tools -NoNewWindow -RedirectStandardOutput (Join-Path $Output 'stdout.log') -RedirectStandardError (Join-Path $Output 'stderr.log') -PassThru
 $nativeHandle=$process.Handle
 $record.test_pid=$process.Id;$record.test_start_time_utc=$process.StartTime.ToUniversalTime().ToString('o')
 $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 while (!$process.HasExited) {
  $process.Refresh();if ($process.HasExited) { break }
  $peak=[Math]::Max($peak,$process.WorkingSet64)
  if ($watch.Elapsed.TotalSeconds -gt 30 -or $process.WorkingSet64 -gt 1073741824) {
   $record.outer_stop=$true;$process.Kill()
   if (!$process.WaitForExit(5000)) { throw 'Owned tool did not exit after guard stop' };break
  }
  Start-Sleep -Milliseconds 50
 }
 if (!$process.WaitForExit(5000)) { throw 'Owned tool did not complete its exit' }
 $record.returncode=[SiriusOwnedExit]::Read($nativeHandle)
 $record.exit_code_source='GetExitCodeProcess on preacquired owned handle after WaitForExit'
} catch { $record.exception=$_.Exception.Message;$record.status='failed' }
finally {
 $record.owned_process_absent=($null -eq $process)
 if ($process) {
  try { if (!$process.HasExited) {$process.Kill();if (!$process.WaitForExit(5000)) {throw 'Owned help tool did not terminate'}} } catch {$errors.Add($_.Exception.Message)}
  try {$record.owned_process_absent=$process.HasExited} catch {$errors.Add($_.Exception.Message)}
 }
 $record.elapsed_seconds=$watch.Elapsed.TotalSeconds;$record.peak_sampled_rss_bytes=$peak
 $record.tool_unchanged=$false
 try {$record.tool_unchanged=((Get-FileHash -LiteralPath $exe -Algorithm SHA256).Hash.ToLowerInvariant() -eq $hash)} catch {$errors.Add('posthash: '+$_.Exception.Message)}
 $record.accessible_remaining_developer_processes=@();$record.developer_census_completed=$false
 try {
  $record.accessible_remaining_developer_processes=@(Get-Process | Where-Object {$_.ProcessName -match '^Radeon(Developer|GPUProfiler)'} | ForEach-Object {@{pid=$_.Id;name=$_.ProcessName;start_utc=$_.StartTime.ToUniversalTime().ToString('o')}})
  $record.developer_census_completed=$true
 } catch {$errors.Add('developer-census: '+$_.Exception.Message)}
 $record.cleanup_errors=$errors.ToArray()
 if ($record.status -ne 'failed') {$record.status=if ($record.outer_stop) {'invalidated'} else {'completed'}}
 if ($errors.Count -or !$record.owned_process_absent -or !$record.tool_unchanged -or !$record.developer_census_completed -or $record.accessible_remaining_developer_processes.Count) {$record.status='failed'}
 try {$record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8} catch {
  $record.status='failed';$record.owner_receipt_error=$_.Exception.Message
  [Console]::Error.WriteLine('Terminal tool owner receipt failed: '+$record.owner_receipt_error)
 }
}
$record | ConvertTo-Json -Depth 8
if ($record.status -ne 'completed' -or $record.returncode -ne 0) {exit 1}
