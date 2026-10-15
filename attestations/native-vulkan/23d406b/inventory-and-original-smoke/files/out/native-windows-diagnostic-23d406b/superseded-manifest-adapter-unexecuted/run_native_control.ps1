param([Parameter(Mandatory=$true)][string]$Stage,[Parameter(Mandatory=$true)][string]$Output,[Parameter(Mandatory=$true)][string]$Filter,[int]$Seconds=320,[long]$RssBytes=4294967296,[string]$Expected=$Filter,[string]$ExecutableName='sirius_backend_tests.exe',[string]$TemporaryDirectory='',[Parameter(Mandatory=$true)][string]$BuildRoot,[Parameter(Mandatory=$true)][string]$ExpectedManifestSha256,[Parameter(Mandatory=$true)][string]$LiveRevision,[string]$RenderDirectory='',[ValidateSet('Control','SystemInventory')][string]$InvocationKind='Control')
$ErrorActionPreference='Stop'
if (!$TemporaryDirectory) { throw 'Task-local temporary directory required' }
if ($TemporaryDirectory) {
 $tempFull=(New-Item -ItemType Directory -Path $TemporaryDirectory -Force).FullName
 $env:TEMP=$tempFull;$env:TMP=$tempFull
 if ($RenderDirectory) { $env:TMPDIR=(New-Item -ItemType Directory -Path $RenderDirectory -Force).FullName } elseif (![string]::IsNullOrEmpty($env:TMPDIR)) { throw 'Unexpected inherited TMPDIR' }
}
$ownerProcess=Get-Process -Id $PID
$ownerHandle=$ownerProcess.Handle
$ownerStart=$ownerProcess.StartTime.ToUniversalTime().ToString('o')
$bootstrap=[ordered]@{schema_version=1;kind='native-owner-bootstrap';owner_pid=$PID;owner_start_time_utc=$ownerStart;temporary_directory=$tempFull;script=$MyInvocation.MyCommand.Path;expected_test_executable=(Join-Path $Stage $ExecutableName);source_revision=$LiveRevision;scope='Early exact birth before Add-Type/preflight; compiler descendant containment remains limited to accessible post-stop observations.'}
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
foreach ($name in @('VK_DRIVER_FILES','VK_ICD_FILENAMES','VK_ADD_DRIVER_FILES','SIRIUS_VULKAN_DEVICE','SIRIUS_PRECISION','SIRIUS_MEMORY_BUDGET_MB','SIRIUS_DISPATCH_TARGET_MS','GTEST_REPEAT','GTEST_TOTAL_SHARDS','GTEST_SHARD_INDEX','GTEST_SHARD_STATUS_FILE','GTEST_FAIL_FAST','GTEST_ALSO_RUN_DISABLED_TESTS')) {
 if (![string]::IsNullOrEmpty([Environment]::GetEnvironmentVariable($name))) { throw ('Unexpected inherited execution override: '+$name) }
}
$busy=@(Get-Process | Where-Object {$_.ProcessName -match '^sirius'})
if ($busy.Count) { throw 'A native Sirius consumer is already live' }
$expectedNames=@($Expected -split ';' | Sort-Object)
if (!$expectedNames.Count -or $Expected -match '[*?]') { throw 'Explicit expected case names required' }
$exe=Join-Path $Stage $ExecutableName
if (!(Test-Path $exe)) { throw 'Missing staged test executable' }
if (Test-Path $Output) { throw 'Output already exists' }
$Output=(New-Item -ItemType Directory -Path $Output).FullName
$Stage=(Get-Item -LiteralPath $Stage).FullName

$exeHash=(Get-FileHash -Algorithm SHA256 -LiteralPath $exe).Hash.ToLowerInvariant()
$BuildRoot=(Get-Item -LiteralPath $BuildRoot).FullName
$manifestPath=Join-Path $BuildRoot 'generated/sirius/windows-diagnostic-build-products.json'
$manifestHash=(Get-FileHash -Algorithm SHA256 -LiteralPath $manifestPath).Hash.ToLowerInvariant()
if ($manifestHash -ne $ExpectedManifestSha256 -or $LiveRevision -cne '23d406b8ca99695f362b88e0c5f9a706f08a7a39') { throw 'Diagnostic product/live source identity mismatch' }
$manifest=Get-Content -Raw -LiteralPath $manifestPath | ConvertFrom-Json
if ($manifest.schema -cne 'sirius-windows-diagnostic-build-products-v1' -or $manifest.source_revision -cne $LiveRevision -or $manifest.repository -cne 'sunyaisblank/Sirius' -or $manifest.run_id -cne '38701994420' -or $manifest.run_attempt -cne '1' -or $manifest.diagnostic_build_only -isnot [bool] -or !$manifest.diagnostic_build_only -or $manifest.qualification_claimed -isnot [bool] -or $manifest.qualification_claimed) { throw 'Wrong diagnostic build-only manifest' }
foreach ($entry in $manifest.files.PSObject.Properties) {
 $key=$entry.Name;$prefix='bin/windows-msvc/'
 if (!$key.StartsWith($prefix,[StringComparison]::Ordinal) -or $key.Contains('\') -or $key.Contains(':')) { throw 'Unsafe diagnostic product key' }
 $relative=$key.Substring($prefix.Length);$parts=@($relative -split '/')
 if (@($parts | Where-Object {$_ -eq '' -or $_ -eq '.' -or $_ -eq '..'}).Count) { throw 'Unsafe diagnostic product path component' }
 $path=[IO.Path]::GetFullPath((Join-Path $BuildRoot $relative))
 if (!$path.StartsWith($BuildRoot.TrimEnd('\')+'\',[StringComparison]::OrdinalIgnoreCase)) { throw 'Diagnostic product escapes build root' }
 if (!(Test-Path -LiteralPath $path -PathType Leaf) -or (Get-Item -LiteralPath $path).Length -ne $entry.Value.bytes -or (Get-FileHash -Algorithm SHA256 -LiteralPath $path).Hash.ToLowerInvariant() -cne $entry.Value.sha256) { throw ('Diagnostic product differs: '+$key) }
}
if ($InvocationKind -eq 'SystemInventory') {
 if ($ExecutableName -cne 'sirius.exe' -or $Filter -cne 'SystemInventory' -or $Expected -cne 'SystemInventory') { throw 'Only exact inventory is admitted' }
 $exeKey='bin/windows-msvc/src/sirius/app/Release/sirius.exe'
} else {
 if ($ExecutableName -cne 'sirius_backend_tests.exe' -or $Filter -cne 'VulkanBackend.SlangKernelMatchesCpuReference' -or $Expected -cne $Filter) { throw 'Only original smoke is admitted' }
 $exeKey='bin/windows-msvc/tests/backend/Release/sirius_backend_tests.exe'
}
$binding=$manifest.files.PSObject.Properties[$exeKey].Value
$expectedPath=[IO.Path]::GetFullPath((Join-Path $BuildRoot $exeKey.Substring('bin/windows-msvc/'.Length)))
if (!$binding -or $exeHash -cne $binding.sha256 -or (Get-Item -LiteralPath $exe).Length -ne $binding.bytes -or ![String]::Equals([IO.Path]::GetFullPath($exe),$expectedPath,[StringComparison]::OrdinalIgnoreCase)) { throw 'Executable differs from exact diagnostic product' }
$loaderPath=Join-Path $Stage 'vulkan-1.dll';$loaderKey=$exeKey.Substring(0,$exeKey.LastIndexOf('/')+1)+'vulkan-1.dll'
$loaderBinding=$manifest.files.PSObject.Properties[$loaderKey].Value
if (!$loaderBinding) { throw 'Exact staged SDK loader is absent from manifest' }
$stdout=Join-Path $Output 'stdout.log';$stderr=Join-Path $Output 'stderr.log';$xml=Join-Path $Output 'gtest.xml'
$arguments=if ($InvocationKind -eq 'SystemInventory') { @('--json','info','system') } else { @("--gtest_filter=$Filter",'--gtest_color=no',('"--gtest_output=xml:'+ $xml+'"')) }
$record=[ordered]@{schema_version=1;kind='bounded-native-Windows-control';start_utc=[DateTime]::UtcNow.ToString('o');owner_pid=$PID;owner_start_time_utc=$ownerStart;test_executable=$exe;test_executable_sha256=$exeHash;arguments=$arguments;working_directory=$Stage;outer_guard_seconds=$Seconds;rss_guard_bytes=$RssBytes;outer_stop=$false;status='running';numerical_control_passed=$false;expected_cases=$expectedNames;qualification_claimed=$false;scope='Exact diagnostic-manifest inventory or original smoke on actual native driver. Diagnostic RSS/time guards retain original scientific checks; no native-build/Mandatory gate, UI, full runtime or release qualification claim.'}
$record.invocation_kind=$InvocationKind
$record.system_inventory_passed=$false
$record.environment=[ordered]@{}
foreach ($name in @('VK_DRIVER_FILES','VK_ICD_FILENAMES','VK_ADD_DRIVER_FILES','SIRIUS_VULKAN_DEVICE','SIRIUS_PRECISION','SIRIUS_MEMORY_BUDGET_MB','SIRIUS_DISPATCH_TARGET_MS','TEMP','TMP','TMPDIR')) { $record.environment[$name]=[Environment]::GetEnvironmentVariable($name) }
$watch=[Diagnostics.Stopwatch]::StartNew();$samples=New-Object 'System.Collections.Generic.List[object]';$process=$null
try {
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
   $loadedLoader=@($bindings | Where-Object {$_.name -ieq 'vulkan-1.dll'})
   if ($loadedLoader.Count -ne 1 -or ![String]::Equals([IO.Path]::GetFullPath($loadedLoader[0].path),[IO.Path]::GetFullPath($loaderPath),[StringComparison]::OrdinalIgnoreCase) -or $loadedLoader[0].sha256 -cne $loaderBinding.sha256 -or $loadedLoader[0].bytes -ne $loaderBinding.bytes) { throw 'Observed loader differs from exact staged SDK loader' }
    ConvertTo-Json -InputObject $bindings -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'driver-modules.json') -Encoding UTF8
   }
  }
  $rss=$process.WorkingSet64
  $samples.Add([pscustomobject]@{elapsed_seconds=$watch.Elapsed.TotalSeconds;working_set_bytes=$rss;private_bytes=$process.PrivateMemorySize64})
  $reason=$null
  if ($rss -gt $RssBytes) { $reason='host_rss' } elseif ($watch.Elapsed.TotalSeconds -gt $Seconds) { $reason='elapsed' }
  if ($reason) { $record.outer_stop=$true;$record.guard_trigger=$reason;$process.Kill();if (!$process.WaitForExit(5000)) { throw 'Guard-stopped owned child did not exit within five seconds' };break }
  if (Test-Path (Join-Path $Output 'driver-modules.json')) { Start-Sleep -Milliseconds 1000 } else { Start-Sleep -Milliseconds 50 }
 }
 if (!$process.WaitForExit(5000)) { throw 'Owned child did not exit within five seconds' };$record.returncode=[SiriusOwnedExit]::Read($nativeHandle);$record.exit_code_source='GetExitCodeProcess on the preacquired owned handle after WaitForExit'
 if ($record.returncode -eq 0 -and !$record.outer_stop -and $InvocationKind -eq 'SystemInventory') {
  $inventory=Get-Content -Raw -LiteralPath $stdout | ConvertFrom-Json
  if (!$inventory.backends.vulkan.available -or $inventory.backends.vulkan.device_count -ne 1 -or $inventory.backends.vulkan.selected_device_name -ne 'AMD Radeon 780M Graphics') { throw 'Unexpected selected native device inventory' }
  $record.system_inventory_passed=$true
 }
 if ($record.returncode -eq 0 -and !$record.outer_stop -and $InvocationKind -eq 'Control') {
  if (!(Test-Path $xml)) { throw 'Completed process did not produce expected GoogleTest XML' }
  [xml]$results=Get-Content -Raw -LiteralPath $xml
  $cases=@($results.SelectNodes('//testcase'))
  $actualNames=@($cases | ForEach-Object { $_.classname+'.'+$_.name } | Sort-Object)
  if (($actualNames -join ';') -ne ($expectedNames -join ';')) { throw 'Executed GoogleTest case names/count differ from expected control' }
  foreach ($suite in $results.SelectNodes('//testsuite')) {
   if ([int]$suite.failures -ne 0 -or [int]$suite.errors -ne 0 -or [int]$suite.skipped -ne 0 -or [int]$suite.disabled -ne 0) { throw 'GoogleTest reported failure/error/skip' }
  }
  if ([int]$results.testsuites.tests -ne $expectedNames.Count -or [int]$results.testsuites.failures -ne 0 -or [int]$results.testsuites.errors -ne 0) { throw 'GoogleTest totals do not establish requested complete control' }
  foreach ($case in $cases) { if ($case.status -ne 'run' -or $case.result -ne 'completed') { throw 'GoogleTest case did not run to completion' } }
  $record.numerical_control_passed=$true
 }
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
 $watch.Stop();$record.end_utc=[DateTime]::UtcNow.ToString('o');$record.elapsed_seconds=$watch.Elapsed.TotalSeconds
 if ($record.status -ne 'failed') { $record.status=if ($record.outer_stop) {'invalidated'} else {'completed'} }
 $record.peak_sampled_rss_bytes=if ($samples.Count) { ($samples.ToArray() | Measure-Object -Property working_set_bytes -Maximum).Maximum } else { 0 }
 $record.source_revision=$manifest.source_revision
 $record.live_source_revision=$LiveRevision
 $record.source_product_manifest_sha256=$manifestHash
 $record.source_product_manifest_path=$manifestPath
 $record.expected_staged_loader=[ordered]@{path=$loaderPath;bytes=$loaderBinding.bytes;sha256=$loaderBinding.sha256}
 $record.native_build_or_mandatory_gate_admitted=$false
 $record.source_product_manifest_unchanged=$false
 try { $record.source_product_manifest_unchanged=((Get-FileHash -Algorithm SHA256 -LiteralPath $manifestPath).Hash.ToLowerInvariant() -eq $manifestHash) } catch { $cleanupErrors.Add('post-exit-gate: '+$_.Exception.Message) }
 if (!$record.source_product_manifest_unchanged) { $cleanupErrors.Add('diagnostic product manifest changed') }
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
if ($record.status -ne 'completed' -or !$record.test_executable_unchanged -or !$record.owned_process_absent -or !$record.source_product_manifest_unchanged -or ($InvocationKind -eq 'Control' -and (!$record.loaded_driver_modules_recorded -or !$record.numerical_control_passed)) -or ($InvocationKind -eq 'SystemInventory' -and !$record.system_inventory_passed)) { exit 1 }
