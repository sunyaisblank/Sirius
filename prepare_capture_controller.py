"""Add owned RDS/RDP capture to the unchanged original numerical control owner."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib
work=Path(__file__).resolve().parent
ps=(work/'run_original_control.ps1').read_text()
assert hashlib.sha256(ps.encode()).hexdigest()=='4394988ece0fdc1cea47681a8aab9bd5e3d4e108f0ad63eb8c678b9219daa760'
def change(old,new):
    global ps
    assert ps.count(old)==1,old
    ps=ps.replace(old,new)
change('[string]$LiveRevision)', '[string]$LiveRevision,[Parameter(Mandatory=$true)][string]$Tools,[Parameter(Mandatory=$true)][string]$PanelHash,[Parameter(Mandatory=$true)][string]$ServiceHash)')
change("$ErrorActionPreference='Stop'", """$ErrorActionPreference='Stop'
$serviceExe=Join-Path $Tools 'RadeonDeveloperServiceCLI.exe'
$panelExe=Join-Path $Tools 'RadeonDeveloperPanelCLI.exe'""")
change("source_revision=$LiveRevision;scope=", "source_revision=$LiveRevision;owned_child_executables=@((Join-Path $Stage $ExecutableName),$serviceExe,$panelExe);scope=")
change("$expectedNames=@($Expected -split ';' | Sort-Object)", """if (@(Get-Process | Where-Object {$_.ProcessName -match '^Radeon(Developer|GPUProfiler)'}).Count) {throw 'Existing developer/profiler consumer'}
if ($Filter -ne 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults' -or $Expected -ne $Filter -or $Seconds -ne 90 -or $RssBytes -ne 4294967296) {throw 'Original one-case scope/guard mismatch'}
foreach ($entry in @(@{exe=$panelExe;hash=$PanelHash},@{exe=$serviceExe;hash=$ServiceHash})) {
 if ((Get-FileHash -Algorithm SHA256 -LiteralPath $entry.exe).Hash.ToLowerInvariant() -ne $entry.hash -or (Get-AuthenticodeSignature -LiteralPath $entry.exe).Status -ne 'Valid') {throw 'Pinned signed profiler tool mismatch'}
}
$expectedNames=@($Expected -split ';' | Sort-Object)""")
change("$stdout=Join-Path $Output 'stdout.log';", "$trace=Join-Path $Output 'original-modern-te.rgp'\n$panelStdout=Join-Path $Output 'panel-stdout.log';$panelStderr=Join-Path $Output 'panel-stderr.log'\n$stdout=Join-Path $Output 'stdout.log';")
change("status='running';numerical_control_passed=", "status='running';scoped_processes=@();capture_client_returncode=$null;capture_trace_sealed=$false;intentional_service_stop=$false;owned_service_absent=$false;owned_capture_client_absent=$false;numerical_control_passed=")
change("Original numerical/reference control on actual native driver. Diagnostic RSS/time guard invalidates incomplete runs; no scientific/device budget changes or full runtime/release qualification claim.", "Original numerical/reference control on actual native driver with dispatch:43:2 instruction/counter trace. Original 90s/4GiB control guard retained. Tool startup/dump and 1GiB sampled tool guards are diagnostic. No clocks, driver experiments, shader instrumentation, performance benefit or full qualification claim.")
change("$watch=[Diagnostics.Stopwatch]::StartNew();$samples=New-Object 'System.Collections.Generic.List[object]';$process=$null", "$watch=[Diagnostics.Stopwatch]::StartNew();$captureWatch=[Diagnostics.Stopwatch]::StartNew();$samples=New-Object 'System.Collections.Generic.List[object]';$process=$null;$service=$null;$panel=$null;$toolPeak=0")
change("try {\n $process=Start-Process", """try {
 $service=Start-Process -FilePath $serviceExe -WorkingDirectory $Tools -NoNewWindow -RedirectStandardOutput (Join-Path $Output 'service-stdout.log') -RedirectStandardError (Join-Path $Output 'service-stderr.log') -PassThru
 $serviceHandle=$service.Handle
 $record.scoped_processes+=@{kind='owned-RDS';pid=$service.Id;birth=$service.StartTime.ToUniversalTime().ToString('o');executable=$serviceExe;sha256=$ServiceHash}
 $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 Start-Sleep -Milliseconds 1000
 if ($service.HasExited) {throw 'Owned service exited before capture client'}
 $panelArguments=@('--mode=profiling','--verbose','--process=sirius_backend_tests.exe','--remote-host=127.0.0.1','--remote-port=27300','--rgp-auto-capture=dispatch:43:2','--rgp-capture-mode=dispatch','--rgp-render-op-count=2','--rgp-instruction-tracing','--rgp-counter-collection','--rgp-sqtt-buffer-size=default',('"--output='+$trace+'"'))
 $record.capture_arguments=$panelArguments
 $panel=Start-Process -FilePath $panelExe -ArgumentList $panelArguments -WorkingDirectory $Tools -NoNewWindow -RedirectStandardOutput $panelStdout -RedirectStandardError $panelStderr -PassThru
 $panelHandle=$panel.Handle
 $record.scoped_processes+=@{kind='owned-RDP-capture';pid=$panel.Id;birth=$panel.StartTime.ToUniversalTime().ToString('o');executable=$panelExe;sha256=$PanelHash}
 $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 $ready=[Diagnostics.Stopwatch]::StartNew();$readyObserved=$false
 while (!$readyObserved) {
  $service.Refresh();$panel.Refresh()
  if ($service.HasExited -or $panel.HasExited) {throw 'Owned capture topology exited before application readiness'}
  $toolPeak=[Math]::Max($toolPeak,[Math]::Max($service.WorkingSet64,$panel.WorkingSet64))
  if ($ready.Elapsed.TotalSeconds -gt 15 -or $toolPeak -gt 1073741824) {throw 'Capture startup readiness/tool RSS bound exceeded; no science launch'}
  $text=Get-Content -Raw -LiteralPath $panelStdout
  $readyObserved=($text -and $text.Contains('Waiting for GPU application to connect...'))
  if (!$readyObserved) {Start-Sleep -Milliseconds 100}
 }
 $record.capture_ready_observed_before_application=$true
 $record.capture_ready_elapsed_seconds=$ready.Elapsed.TotalSeconds
 $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 $watch.Restart()
 $process=Start-Process""")
change("  $rss=$process.WorkingSet64", """  $service.Refresh();$panel.Refresh()
  if ($service.HasExited) {throw 'Owned service exited before capture completion'}
  if ($panel.HasExited) {
   if ([SiriusOwnedExit]::Read($panelHandle) -ne 0) {throw 'Capture client failed during original control'}
  } else {
   $toolPeak=[Math]::Max($toolPeak,$panel.WorkingSet64)
  }
  $toolPeak=[Math]::Max($toolPeak,$service.WorkingSet64)
  if ($toolPeak -gt 1073741824) {throw 'Diagnostic profiler tool RSS bound exceeded'}
  $rss=$process.WorkingSet64""")
change(" } catch {\n $record.status='failed'", """ $record.science_elapsed_seconds=$watch.Elapsed.TotalSeconds
 if (!$record.numerical_control_passed) {throw 'Incomplete original numerical control cannot establish a trace result'}
 $dump=[Diagnostics.Stopwatch]::StartNew()
 while (!$panel.HasExited) {
  $panel.Refresh();$service.Refresh()
  if ($service.HasExited -or $dump.Elapsed.TotalSeconds -gt 30) {throw 'Capture dump/service lifetime bound exceeded'}
  $toolPeak=[Math]::Max($toolPeak,$service.WorkingSet64)
  if (!$panel.HasExited) {$toolPeak=[Math]::Max($toolPeak,$panel.WorkingSet64)}
  if ($toolPeak -gt 1073741824) {throw 'Diagnostic capture dump RSS bound exceeded'}
  Start-Sleep -Milliseconds 100
 }
 if (!$panel.WaitForExit(5000)) {throw 'Owned capture client did not complete exit'}
 $record.capture_client_returncode=[SiriusOwnedExit]::Read($panelHandle)
 if ($record.capture_client_returncode -ne 0) {throw 'Capture client did not complete successfully'}
 $size=(Get-Item -LiteralPath $trace).Length
 if ($size -le 0 -or $size -gt 268435456) {throw 'Trace missing or exceeds bounded analysis input size'}
 $record.trace=@{path=$trace;bytes=$size;sha256=(Get-FileHash -LiteralPath $trace -Algorithm SHA256).Hash.ToLowerInvariant()}
 $record.capture_trace_sealed=$true
 } catch {
 $record.status='failed'""")
change(" $watch.Stop();$record.end_utc=", """ $record.owned_capture_client_absent=($null -eq $panel)
 if ($panel) {
  try {
   if (!$panel.HasExited) {$record.intentional_failed_capture_stop=$true;$panel.Kill();if (!$panel.WaitForExit(5000)) {throw 'Owned capture client did not terminate'}}
   $record.capture_client_terminal_returncode=[SiriusOwnedExit]::Read($panelHandle)
   $record.owned_capture_client_absent=$panel.HasExited
  } catch {$cleanupErrors.Add('capture-client-stop: '+$_.Exception.Message)}
 }
 $record.owned_service_absent=($null -eq $service)
 if ($service) {
  try {
   if ($service.HasExited) {$cleanupErrors.Add('Owned service exited before deliberate stop')} else {$record.intentional_service_stop=$true;$service.Kill();if (!$service.WaitForExit(5000)) {throw 'Owned service did not terminate'}}
   $record.service_terminal_returncode=[SiriusOwnedExit]::Read($serviceHandle)
   $record.owned_service_absent=$service.HasExited
  } catch {$cleanupErrors.Add('service-stop: '+$_.Exception.Message)}
 }
 $record.profiler_tools_unchanged=$false;$record.developer_census_completed=$false;$record.accessible_remaining_developer_processes=@()
 try {$record.profiler_tools_unchanged=((Get-FileHash -LiteralPath $panelExe -Algorithm SHA256).Hash.ToLowerInvariant() -eq $PanelHash -and (Get-FileHash -LiteralPath $serviceExe -Algorithm SHA256).Hash.ToLowerInvariant() -eq $ServiceHash)} catch {$cleanupErrors.Add('profiler-posthash: '+$_.Exception.Message)}
 try {$record.accessible_remaining_developer_processes=@(Get-Process | Where-Object {$_.ProcessName -match '^Radeon(Developer|GPUProfiler)'} | ForEach-Object {@{pid=$_.Id;name=$_.ProcessName;start_utc=$_.StartTime.ToUniversalTime().ToString('o')}});$record.developer_census_completed=$true} catch {$cleanupErrors.Add('profiler-census: '+$_.Exception.Message)}
 if (!$record.owned_capture_client_absent -or !$record.owned_service_absent -or !$record.profiler_tools_unchanged -or !$record.developer_census_completed -or $record.accessible_remaining_developer_processes.Count) {$cleanupErrors.Add('Profiler closure/identity/census incomplete')}
 $record.sampled_peak_profiler_process_rss_bytes=$toolPeak;$record.capture_owner_elapsed_seconds=$captureWatch.Elapsed.TotalSeconds
 $watch.Stop();$record.end_utc=""")
change("if ($record.status -ne 'completed' -or !$record.loaded_driver_modules_recorded", "if ($record.status -ne 'completed' -or !$record.capture_trace_sealed -or $record.capture_client_returncode -ne 0 -or !$record.capture_ready_observed_before_application -or !$record.intentional_service_stop -or !$record.owned_service_absent -or !$record.owned_capture_client_absent -or !$record.loaded_driver_modules_recorded")
with (work/'run_capture_control.ps1').open('x') as f:f.write(ps)
print(len(ps.encode()),hashlib.sha256(ps.encode()).hexdigest())
