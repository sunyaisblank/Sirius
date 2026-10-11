param([Parameter(Mandatory=$true)][string]$Tools,[Parameter(Mandatory=$true)][string]$Output,[Parameter(Mandatory=$true)][ValidateSet('service-info-loopback')][string]$Action,[Parameter(Mandatory=$true)][string]$ExpectedHash,[Parameter(Mandatory=$true)][string]$ExpectedServiceHash,[Parameter(Mandatory=$true)][string]$LiveRevision)
$ErrorActionPreference='Stop'
if (Test-Path -LiteralPath $Output) { throw 'One-shot tool output already exists' }
$Output=(New-Item -ItemType Directory -Path $Output).FullName
$tempFull=(New-Item -ItemType Directory -Path (Join-Path $Output 'temp')).FullName
$env:TEMP=$tempFull;$env:TMP=$tempFull
$tool='RadeonDeveloperPanelCLI.exe'
$argument=@('--system-info','--verbose','--process','sirius_backend_tests.exe','--remote-host','127.0.0.1','--remote-port','27300')
$exe=Join-Path $Tools $tool
$serviceExe=Join-Path $Tools 'RadeonDeveloperServiceCLI.exe'
$ownerProcess=Get-Process -Id $PID;$ownerHandle=$ownerProcess.Handle
$ownerStart=$ownerProcess.StartTime.ToUniversalTime().ToString('o')
$bootstrap=[ordered]@{owner_pid=$PID;owner_start_time_utc=$ownerStart;expected_test_executable=$exe;owned_child_executables=@($exe,$serviceExe);scope='Early exact service/client owner birth. Read-only system information only; no scientific workload or capture.'}
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
if ($existing.Count -or @(Get-Process | Where-Object {$_.ProcessName -match '^sirius'}).Count) { throw 'Existing developer/profiler/Sirius consumer: service diagnostic declined' }
$hash=(Get-FileHash -LiteralPath $exe -Algorithm SHA256).Hash.ToLowerInvariant()
if ($hash -ne $ExpectedHash -or $LiveRevision -notmatch '^[0-9a-f]{40}$') { throw 'Pinned tool/source mismatch' }
$signature=Get-AuthenticodeSignature -LiteralPath $exe
$sig=[ordered]@{status=$signature.Status.ToString();subject=$signature.SignerCertificate.Subject;thumbprint=$signature.SignerCertificate.Thumbprint;file_version=(Get-Item -LiteralPath $exe).VersionInfo.FileVersion;product_version=(Get-Item -LiteralPath $exe).VersionInfo.ProductVersion}
$sig | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'signature.json') -Encoding UTF8
if ($signature.Status -ne 'Valid') { throw 'Official tool signature did not verify' }
$serviceHash=(Get-FileHash -LiteralPath $serviceExe -Algorithm SHA256).Hash.ToLowerInvariant()
if ($serviceHash -ne $ExpectedServiceHash) {throw 'Pinned service mismatch'}
$serviceSignature=Get-AuthenticodeSignature -LiteralPath $serviceExe
@{status=$serviceSignature.Status.ToString();subject=$serviceSignature.SignerCertificate.Subject;thumbprint=$serviceSignature.SignerCertificate.Thumbprint;file_version=(Get-Item -LiteralPath $serviceExe).VersionInfo.FileVersion} | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'service-signature.json') -Encoding UTF8
if ($serviceSignature.Status -ne 'Valid') {throw 'Official service signature did not verify'}
$record=[ordered]@{action=$Action;owner_pid=$PID;owner_start_time_utc=$ownerStart;test_executable=$exe;test_executable_sha256=$hash;arguments=@($argument);live_source_revision=$LiveRevision;outer_guard_seconds=30;rss_guard_bytes=1073741824;outer_stop=$false;status='running';scoped_processes=@();owned_service_absent=$false;intentional_service_stop=$false;scope='Pinned owned service and system-info client only. Deliberate owned service stop after information request; no application capture, clock/driver experiment or scientific qualification.'}
$process=$null;$service=$null;$watch=[Diagnostics.Stopwatch]::StartNew();$peak=0;$errors=New-Object 'System.Collections.Generic.List[string]'
try {
 $service=Start-Process -FilePath $serviceExe -WorkingDirectory $Tools -NoNewWindow -RedirectStandardOutput (Join-Path $Output 'service-stdout.log') -RedirectStandardError (Join-Path $Output 'service-stderr.log') -PassThru
 $serviceHandle=$service.Handle
 $record.scoped_processes=@(@{kind='owned-RDS';pid=$service.Id;birth=$service.StartTime.ToUniversalTime().ToString('o');executable=$serviceExe;sha256=$serviceHash})
 $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 Start-Sleep -Milliseconds 1000
 if ($service.HasExited) {throw 'Owned service exited before the information client'}
 $process=Start-Process -FilePath $exe -ArgumentList $argument -WorkingDirectory $Tools -NoNewWindow -RedirectStandardOutput (Join-Path $Output 'stdout.log') -RedirectStandardError (Join-Path $Output 'stderr.log') -PassThru
 $nativeHandle=$process.Handle
 $record.test_pid=$process.Id;$record.test_start_time_utc=$process.StartTime.ToUniversalTime().ToString('o')
 $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 while (!$process.HasExited) {
  $process.Refresh();if ($process.HasExited) { break }
  $service.Refresh();if ($service.HasExited) {throw 'Owned service exited during information request'}
  $peak=[Math]::Max($peak,[Math]::Max($process.WorkingSet64,$service.WorkingSet64))
  if ($watch.Elapsed.TotalSeconds -gt 30 -or $process.WorkingSet64 -gt 1073741824 -or $service.WorkingSet64 -gt 1073741824) {
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
 $record.owned_service_absent=($null -eq $service)
 if ($service) {
  try {
   if (!$service.HasExited) {$record.intentional_service_stop=$true;$service.Kill();if (!$service.WaitForExit(5000)) {throw 'Owned service did not terminate'}} else {$errors.Add('Owned service exited before its deliberate stop')}
   $record.service_returncode=[SiriusOwnedExit]::Read($serviceHandle)
   $record.owned_service_absent=$service.HasExited
  } catch {$errors.Add('service-stop: '+$_.Exception.Message)}
 }
 $record.service_unchanged=$false
 try {$record.service_unchanged=((Get-FileHash -LiteralPath $serviceExe -Algorithm SHA256).Hash.ToLowerInvariant() -eq $serviceHash)} catch {$errors.Add('service-posthash: '+$_.Exception.Message)}
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
 if (!$record.owned_service_absent -or !$record.service_unchanged -or $errors.Count -or !$record.owned_process_absent -or !$record.tool_unchanged -or !$record.developer_census_completed -or $record.accessible_remaining_developer_processes.Count) {$record.status='failed'}
 try {$record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8} catch {
  $record.status='failed';$record.owner_receipt_error=$_.Exception.Message
  [Console]::Error.WriteLine('Terminal tool owner receipt failed: '+$record.owner_receipt_error)
 }
}
$record | ConvertTo-Json -Depth 8
if ($record.status -ne 'completed' -or $record.returncode -ne 0) {exit 1}
