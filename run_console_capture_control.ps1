param([Parameter(Mandatory=$true)][string]$Stage,[Parameter(Mandatory=$true)][string]$Output,[Parameter(Mandatory=$true)][string]$Filter,[int]$Seconds=320,[long]$RssBytes=4294967296,[string]$Expected=$Filter,[string]$ExecutableName='sirius_backend_tests.exe',[string]$TemporaryDirectory='',[Parameter(Mandatory=$true)][string]$BuildRoot,[Parameter(Mandatory=$true)][string]$ExpectedGateSha256,[Parameter(Mandatory=$true)][string]$LiveRevision,[Parameter(Mandatory=$true)][string]$Tools,[Parameter(Mandatory=$true)][string]$PanelHash,[Parameter(Mandatory=$true)][string]$ServiceHash,[Parameter(Mandatory=$true)][int]$RequiredConsoleWidth)
$ErrorActionPreference='Stop'
$serviceExe=Join-Path $Tools 'RadeonDeveloperServiceCLI.exe'
$panelExe=Join-Path $Tools 'RadeonDeveloperPanelCLI.exe'
if (!$TemporaryDirectory) { throw 'Task-local temporary directory required' }
if ($TemporaryDirectory) {
 $tempFull=(New-Item -ItemType Directory -Path $TemporaryDirectory -Force).FullName
 $env:TEMP=$tempFull;$env:TMP=$tempFull
}
$ownerProcess=Get-Process -Id $PID
$ownerHandle=$ownerProcess.Handle
$ownerStart=$ownerProcess.StartTime.ToUniversalTime().ToString('o')
$bootstrap=[ordered]@{schema_version=1;kind='native-owner-bootstrap';owner_pid=$PID;owner_start_time_utc=$ownerStart;temporary_directory=$tempFull;script=$MyInvocation.MyCommand.Path;expected_test_executable=(Join-Path $Stage $ExecutableName);source_revision=$LiveRevision;owned_child_executables=@((Join-Path $Stage $ExecutableName),$serviceExe,$panelExe);scope='Early exact birth before Add-Type/preflight; compiler descendant containment remains limited to accessible post-stop observations.'}
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
Add-Type -TypeDefinition (Get-Content -Raw -LiteralPath (Join-Path $PSScriptRoot 'owned_console_capture.cs'))
foreach ($name in @('VK_DRIVER_FILES','VK_ICD_FILENAMES','VK_ADD_DRIVER_FILES','SIRIUS_VULKAN_DEVICE','SIRIUS_PRECISION','SIRIUS_MEMORY_BUDGET_MB','SIRIUS_DISPATCH_TARGET_MS','GTEST_REPEAT','GTEST_TOTAL_SHARDS','GTEST_SHARD_INDEX','GTEST_SHARD_STATUS_FILE','GTEST_FAIL_FAST','GTEST_ALSO_RUN_DISABLED_TESTS')) {
 if (![string]::IsNullOrEmpty([Environment]::GetEnvironmentVariable($name))) { throw ('Unexpected inherited execution override: '+$name) }
}
$busy=@(Get-Process | Where-Object {$_.ProcessName -match '^sirius'})
if ($busy.Count) { throw 'A native Sirius consumer is already live' }
if (@(Get-Process | Where-Object {$_.ProcessName -match '^Radeon(Developer|GPUProfiler)'}).Count) {throw 'Existing developer/profiler consumer'}
if ($Filter -ne 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults' -or $Expected -ne $Filter -or $Seconds -ne 90 -or $RssBytes -ne 4294967296) {throw 'Original one-case scope/guard mismatch'}
foreach ($entry in @(@{exe=$panelExe;hash=$PanelHash},@{exe=$serviceExe;hash=$ServiceHash})) {
 if ((Get-FileHash -Algorithm SHA256 -LiteralPath $entry.exe).Hash.ToLowerInvariant() -ne $entry.hash -or (Get-AuthenticodeSignature -LiteralPath $entry.exe).Status -ne 'Valid') {throw 'Pinned signed profiler tool mismatch'}
}
$expectedNames=@($Expected -split ';' | Sort-Object)
if (!$expectedNames.Count -or $Expected -match '[*?]') { throw 'Explicit expected case names required' }
$exe=Join-Path $Stage $ExecutableName
if (!(Test-Path $exe)) { throw 'Missing staged test executable' }
if (Test-Path $Output) { throw 'Output already exists' }
$Output=(New-Item -ItemType Directory -Path $Output).FullName
$Stage=(Get-Item -LiteralPath $Stage).FullName

$exeHash=(Get-FileHash -Algorithm SHA256 -LiteralPath $exe).Hash.ToLowerInvariant()
$BuildRoot=(Get-Item -LiteralPath $BuildRoot).FullName
$gatePath=Join-Path $BuildRoot 'generated/sirius/native_build_gate.json'
$gateHash=(Get-FileHash -Algorithm SHA256 -LiteralPath $gatePath).Hash.ToLowerInvariant()
if ($gateHash -ne $ExpectedGateSha256 -or $LiveRevision -notmatch '^[0-9a-f]{40}$') { throw 'Native producer/live source gate identity mismatch' }
$gate=Get-Content -Raw -LiteralPath $gatePath | ConvertFrom-Json
if (!$gate.source.clean -or $gate.source.revision -notmatch '^[0-9a-f]{40}$') { throw 'Invalid native producer source receipt' }
$logicalName=[IO.Path]::GetFileNameWithoutExtension($ExecutableName)
$binding=$gate.tested_artifacts.$logicalName
$expectedPath=[IO.Path]::GetFullPath((Join-Path $BuildRoot $binding.path))
if ($binding.root -ne 'build' -or $exeHash -ne $binding.sha256 -or (Get-Item -LiteralPath $exe).Length -ne $binding.bytes -or [IO.Path]::GetFullPath($exe) -ne $expectedPath) { throw 'Executable does not match exact native producer gate' }
$trace=Join-Path $Output 'original-modern-te.rgp'
$panelScreenText=Join-Path $Output 'panel-screen.txt'
$stdout=Join-Path $Output 'stdout.log';$stderr=Join-Path $Output 'stderr.log';$xml=Join-Path $Output 'gtest.xml'
$arguments=@("--gtest_filter=$Filter",'--gtest_color=no',('"--gtest_output=xml:'+ $xml+'"'))
$record=[ordered]@{schema_version=1;kind='bounded-native-Windows-control';start_utc=[DateTime]::UtcNow.ToString('o');owner_pid=$PID;owner_start_time_utc=$ownerStart;test_executable=$exe;test_executable_sha256=$exeHash;arguments=$arguments;working_directory=$Stage;outer_guard_seconds=$Seconds;rss_guard_bytes=$RssBytes;outer_stop=$false;status='running';scoped_processes=@();capture_client_returncode=$null;capture_trace_sealed=$false;intentional_service_stop=$false;owned_service_absent=$false;owned_capture_client_absent=$false;numerical_control_passed=$false;expected_cases=$expectedNames;qualification_claimed=$false;scope='Original numerical/reference control on actual native driver with dispatch:43:2 instruction/counter trace. Original 90s/4GiB control guard retained. Tool startup/dump and 1GiB sampled tool guards are diagnostic. No clocks, driver experiments, shader instrumentation, performance benefit or full qualification claim.'}
$record.environment=[ordered]@{}
foreach ($name in @('VK_DRIVER_FILES','VK_ICD_FILENAMES','VK_ADD_DRIVER_FILES','SIRIUS_VULKAN_DEVICE','SIRIUS_PRECISION','SIRIUS_MEMORY_BUDGET_MB','SIRIUS_DISPATCH_TARGET_MS','TEMP','TMP')) { $record.environment[$name]=[Environment]::GetEnvironmentVariable($name) }
$watch=[Diagnostics.Stopwatch]::StartNew();$captureWatch=[Diagnostics.Stopwatch]::StartNew();$samples=New-Object 'System.Collections.Generic.List[object]';$process=$null;$service=$null;$panel=$null;$toolPeak=0;$consoleAttached=$false;$lastScreen=$null
try {
 $service=Start-Process -FilePath $serviceExe -WorkingDirectory $Tools -NoNewWindow -RedirectStandardOutput (Join-Path $Output 'service-stdout.log') -RedirectStandardError (Join-Path $Output 'service-stderr.log') -PassThru
 $serviceHandle=$service.Handle
 $record.scoped_processes+=@{kind='owned-RDS';pid=$service.Id;birth=$service.StartTime.ToUniversalTime().ToString('o');executable=$serviceExe;sha256=$ServiceHash}
 $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 Start-Sleep -Milliseconds 1000
 if ($service.HasExited) {throw 'Owned service exited before capture client'}
 $panelArguments=@('--mode=profiling','--verbose','--process=sirius_backend_tests.exe','--remote-host=127.0.0.1','--remote-port=27300','--rgp-auto-capture=dispatch:43:2','--rgp-capture-mode=dispatch','--rgp-render-op-count=2','--rgp-instruction-tracing','--rgp-counter-collection','--rgp-sqtt-buffer-size=default',('"--output='+$trace+'"'))
 $record.capture_arguments=$panelArguments
 [SiriusOwnedConsole]::ConfigureWidth($RequiredConsoleWidth)
 $record.required_console_width=$RequiredConsoleWidth
 $record.capture_output_kind='Counted owned-console screen cells; combined channels, not stdout byte history'
 $consolePid=[SiriusOwnedConsole]::StartSuspended($panelExe,($panelArguments -join ' '),$Tools)
 $panel=Get-Process -Id $consolePid
 $panelHandle=$panel.Handle
 $record.scoped_processes+=@{kind='owned-RDP-capture';pid=$panel.Id;birth=$panel.StartTime.ToUniversalTime().ToString('o');executable=$panelExe;sha256=$PanelHash}
 $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 [SiriusOwnedConsole]::Resume()
 $attach=[Diagnostics.Stopwatch]::StartNew()
 while (!$consoleAttached) {
  $service.Refresh();$panel.Refresh()
  $toolPeak=[Math]::Max($toolPeak,[Math]::Max($service.WorkingSet64,$panel.WorkingSet64))
  if ($service.HasExited -or $panel.HasExited -or $attach.Elapsed.TotalSeconds -gt 3 -or $toolPeak -gt 1073741824) {throw 'Owned running-console attachment/lifetime/RSS bound exceeded; no science launch'}
  $consoleAttached=[SiriusOwnedConsole]::TryAttach([uint32]$panel.Id)
  $record.last_console_attach_error=[SiriusOwnedConsole]::LastAttachError
  if (!$consoleAttached) {Start-Sleep -Milliseconds 50}
 }
 $lastScreen=[SiriusOwnedConsole]::Read([uint32]$PID,[uint32]$panel.Id,0)
 $lastScreen | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'initial-screen.json') -Encoding UTF8
 $record.console_geometry_verified_after_attachment=$true
 $ready=[Diagnostics.Stopwatch]::StartNew();$readyObserved=$false
 while (!$readyObserved) {
  $service.Refresh();$panel.Refresh()
  if ($service.HasExited -or $panel.HasExited) {throw 'Owned capture topology exited before application readiness'}
  $toolPeak=[Math]::Max($toolPeak,[Math]::Max($service.WorkingSet64,$panel.WorkingSet64))
  if ($ready.Elapsed.TotalSeconds -gt 15 -or $toolPeak -gt 1073741824) {throw 'Capture startup readiness/tool RSS bound exceeded; no science launch'}
  $lastScreen=[SiriusOwnedConsole]::Read([uint32]$PID,[uint32]$panel.Id,0)
  $readyObserved=($lastScreen.text.Contains('Profiling (RGP) feature enabled (auto-capture enabled)') -and $lastScreen.text.Contains('Waiting for GPU application to connect...'))
  if (!$readyObserved) {Start-Sleep -Milliseconds 100}
 }
 $lastScreen | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'ready-screen.json') -Encoding UTF8
 $record.capture_ready_observed_before_application=$true
 $record.capture_ready_elapsed_seconds=$ready.Elapsed.TotalSeconds
 $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 $watch.Restart()
 $process=Start-Process -FilePath $exe -ArgumentList $arguments -WorkingDirectory $Stage -NoNewWindow -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
 $nativeHandle=$process.Handle
 $record.test_pid=$process.Id;$record.test_start_time_utc=$process.StartTime.ToUniversalTime().ToString('o')
 $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 while (!$process.HasExited) {
  $process.Refresh()
  if ($process.HasExited) { break }
  if (!(Test-Path (Join-Path $Output 'driver-modules.json'))) {
   $modules=@($process.Modules | Where-Object {$_.ModuleName -match '(?i)amdvlk|^vulkan-1\.dll$'})
   if (@($modules | Where-Object {$_.ModuleName -match '(?i)amdvlk'}).Count -gt 0 -and @($modules | Where-Object {$_.ModuleName -ieq 'vulkan-1.dll'}).Count -gt 0) {
    $bindings=@($modules | ForEach-Object {[ordered]@{name=$_.ModuleName;path=$_.FileName;bytes=(Get-Item -LiteralPath $_.FileName).Length;sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $_.FileName).Hash.ToLowerInvariant();version=$_.FileVersionInfo.FileVersion}})
    ConvertTo-Json -InputObject $bindings -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'driver-modules.json') -Encoding UTF8
   }
  }
  $service.Refresh();$panel.Refresh()
  if ($service.HasExited) {throw 'Owned service exited before capture completion'}
  if ($panel.HasExited) {
   if ([SiriusOwnedExit]::Read($panelHandle) -ne 0) {throw 'Capture client failed during original control'}
  } else {
   $toolPeak=[Math]::Max($toolPeak,$panel.WorkingSet64)
  }
  $toolPeak=[Math]::Max($toolPeak,$service.WorkingSet64)
  if ($toolPeak -gt 1073741824) {throw 'Diagnostic profiler tool RSS bound exceeded'}
  $rss=$process.WorkingSet64
  $samples.Add([pscustomobject]@{elapsed_seconds=$watch.Elapsed.TotalSeconds;working_set_bytes=$rss;private_bytes=$process.PrivateMemorySize64})
  $reason=$null
  if ($rss -gt $RssBytes) { $reason='host_rss' } elseif ($watch.Elapsed.TotalSeconds -gt $Seconds) { $reason='elapsed' }
  if ($reason) { $record.outer_stop=$true;$record.guard_trigger=$reason;$process.Kill();if (!$process.WaitForExit(5000)) { throw 'Guard-stopped owned child did not exit within five seconds' };break }
  if (Test-Path (Join-Path $Output 'driver-modules.json')) { Start-Sleep -Milliseconds 1000 } else { Start-Sleep -Milliseconds 50 }
 }
 if (!$process.WaitForExit(5000)) { throw 'Owned child did not exit within five seconds' };$record.returncode=[SiriusOwnedExit]::Read($nativeHandle);$record.exit_code_source='GetExitCodeProcess on the preacquired owned handle after WaitForExit'
 if ($record.returncode -eq 0 -and !$record.outer_stop) {
  if (!(Test-Path $xml)) { throw 'Completed process did not produce expected GoogleTest XML' }
  [xml]$results=Get-Content -Raw -LiteralPath $xml
  $cases=@($results.SelectNodes('//testcase'))
  $actualNames=@($cases | ForEach-Object { $_.classname+'.'+$_.name } | Sort-Object)
  if (($actualNames -join ';') -ne ($expectedNames -join ';')) { throw 'Executed GoogleTest case names/count differ from expected control' }
  foreach ($suite in $results.SelectNodes('//testsuite')) {
   if ([int]$suite.failures -ne 0 -or [int]$suite.errors -ne 0 -or [int]$suite.skipped -ne 0) { throw 'GoogleTest reported failure/error/skip' }
  }
  if ([int]$results.testsuites.tests -ne $expectedNames.Count -or [int]$results.testsuites.failures -ne 0 -or [int]$results.testsuites.errors -ne 0) { throw 'GoogleTest totals do not establish requested complete control' }
  $record.numerical_control_passed=$true
 }
 $record.science_elapsed_seconds=$watch.Elapsed.TotalSeconds
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
 $lastScreen=[SiriusOwnedConsole]::Read([uint32]$PID,[uint32]$panel.Id,[uint32]$process.Id)
 $lastScreen | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'final-screen.json') -Encoding UTF8
 $lastScreen.text | Set-Content -LiteralPath $panelScreenText -Encoding UTF8
 $record.capture_client_returncode=[SiriusOwnedExit]::Read($panelHandle)
 if ($record.capture_client_returncode -ne 0) {throw 'Capture client did not complete successfully'}
 $size=(Get-Item -LiteralPath $trace).Length
 if ($size -le 0 -or $size -gt 268435456) {throw 'Trace missing or exceeds bounded analysis input size'}
 $record.trace=@{path=$trace;bytes=$size;sha256=(Get-FileHash -LiteralPath $trace -Algorithm SHA256).Hash.ToLowerInvariant()}
 $record.capture_trace_sealed=$true
 } catch {
 $record.status='failed';$record.exception=$_.Exception.Message
} finally {
 $cleanupErrors=New-Object 'System.Collections.Generic.List[string]'
 $record.owned_process_absent=($null -eq $process)
 if ($process) {
  try {
   if (!$process.HasExited) {
    $process.Kill()
    if (!$process.WaitForExit(5000)) { throw 'Owned child did not terminate after kill' }
   }
  } catch { $cleanupErrors.Add('termination: '+$_.Exception.Message) }
  try { $record.owned_process_absent=$process.HasExited } catch { $cleanupErrors.Add('terminal-state: '+$_.Exception.Message) }
 }
 $record.owned_capture_client_absent=($null -eq $panel)
 if ($panel) {
  try {
   if (!$panel.HasExited) {$record.intentional_failed_capture_stop=$true;$panel.Kill();if (!$panel.WaitForExit(5000)) {throw 'Owned capture client did not terminate'}}
   $record.capture_client_terminal_returncode=[SiriusOwnedExit]::Read($panelHandle)
   $record.owned_capture_client_absent=$panel.HasExited
  } catch {$cleanupErrors.Add('capture-client-stop: '+$_.Exception.Message)}
 }
 $record.console_detached=$false
 try {if ($lastScreen) {$lastScreen | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'last-screen.json') -Encoding UTF8}} catch {$cleanupErrors.Add('last-screen-receipt: '+$_.Exception.Message)}
 try {[SiriusOwnedConsole]::Detach();$record.console_detached=$true} catch {$cleanupErrors.Add('console-detach: '+$_.Exception.Message)}
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
 if (!$record.console_detached -or !$record.owned_capture_client_absent -or !$record.owned_service_absent -or !$record.profiler_tools_unchanged -or !$record.developer_census_completed -or $record.accessible_remaining_developer_processes.Count) {$cleanupErrors.Add('Profiler closure/identity/census incomplete')}
 $record.sampled_peak_profiler_process_rss_bytes=$toolPeak;$record.capture_owner_elapsed_seconds=$captureWatch.Elapsed.TotalSeconds
 $watch.Stop();$record.end_utc=[DateTime]::UtcNow.ToString('o');$record.elapsed_seconds=$watch.Elapsed.TotalSeconds
 if ($record.status -ne 'failed') { $record.status=if ($record.outer_stop) {'invalidated'} else {'completed'} }
 $record.peak_sampled_rss_bytes=if ($samples.Count) { ($samples.ToArray() | Measure-Object -Property working_set_bytes -Maximum).Maximum } else { 0 }
 $record.source_revision=$gate.source.revision
 $record.live_source_revision=$LiveRevision
 $record.source_build_gate_sha256=$gateHash
 $record.source_build_gate_unchanged=$false
 try { $record.source_build_gate_unchanged=((Get-FileHash -Algorithm SHA256 -LiteralPath $gatePath).Hash.ToLowerInvariant() -eq $gateHash) } catch { $cleanupErrors.Add('post-exit-gate: '+$_.Exception.Message) }
 if (!$record.source_build_gate_unchanged) { $cleanupErrors.Add('native producer gate changed') }
 $record.loaded_driver_modules_recorded=$false
 try { $record.loaded_driver_modules_recorded=(Test-Path (Join-Path $Output 'driver-modules.json')) } catch { $cleanupErrors.Add('module-receipt: '+$_.Exception.Message) }
 $record.test_executable_unchanged=$false
 try { $record.test_executable_unchanged=((Get-FileHash -Algorithm SHA256 -LiteralPath $exe).Hash.ToLowerInvariant() -eq $exeHash) } catch { $cleanupErrors.Add('post-exit-hash: '+$_.Exception.Message) }
 try { $samples.ToArray() | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'memory.json') -Encoding UTF8 } catch { $cleanupErrors.Add('memory-receipt: '+$_.Exception.Message) }
 if ($cleanupErrors.Count -or !$record.owned_process_absent) { $record.status='failed' }
 $record.cleanup_errors=$cleanupErrors.ToArray()
 try { $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8 } catch {
  $record.status='failed';$record.owner_receipt_error=$_.Exception.Message
  [Console]::Error.WriteLine('Terminal owner receipt failed: '+$record.owner_receipt_error)
 }
}
[pscustomobject]@{status=$record.status;returncode=$record.returncode;elapsed_seconds=$record.elapsed_seconds;peak_sampled_rss_bytes=$record.peak_sampled_rss_bytes;owned_process_absent=$record.owned_process_absent} | ConvertTo-Json -Compress
if ($record.status -ne 'completed' -or !$record.console_detached -or !$record.capture_trace_sealed -or $record.capture_client_returncode -ne 0 -or !$record.capture_ready_observed_before_application -or !$record.intentional_service_stop -or !$record.owned_service_absent -or !$record.owned_capture_client_absent -or !$record.loaded_driver_modules_recorded -or !$record.numerical_control_passed -or !$record.test_executable_unchanged -or !$record.owned_process_absent) { exit 1 }
