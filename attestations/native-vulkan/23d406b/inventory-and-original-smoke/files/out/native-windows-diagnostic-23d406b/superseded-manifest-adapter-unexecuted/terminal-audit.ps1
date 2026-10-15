param([Parameter(Mandatory=$true)][string]$Workspace,[string]$RunName='native-run',[string]$Receipt='native-terminal-audit.json',[Parameter(Mandatory=$true)][string]$BuildRoot)
$ErrorActionPreference='Stop'
$owner=Get-Content -Raw -LiteralPath (Join-Path $Workspace ($RunName+'\owner.json')) | ConvertFrom-Json
$manifestPath=Join-Path $BuildRoot 'generated/sirius/windows-diagnostic-build-products.json'
$manifestHash=(Get-FileHash -Algorithm SHA256 -LiteralPath $manifestPath).Hash.ToLowerInvariant()
if ($manifestHash -cne $owner.source_product_manifest_sha256 -or !$owner.source_product_manifest_unchanged -or $owner.native_build_or_mandatory_gate_admitted) { throw 'Diagnostic product authority changed' }
$manifest=Get-Content -Raw -LiteralPath $manifestPath | ConvertFrom-Json
if ($manifest.schema -cne 'sirius-windows-diagnostic-build-products-v1' -or $manifest.source_revision -cne '23d406b8ca99695f362b88e0c5f9a706f08a7a39' -or $manifest.repository -cne 'sunyaisblank/Sirius' -or $manifest.run_id -cne '38701994420' -or $manifest.run_attempt -cne '1' -or !$manifest.diagnostic_build_only -or $manifest.qualification_claimed) { throw 'Wrong diagnostic manifest' }
$exeKey=if ($owner.invocation_kind -eq 'SystemInventory') {'bin/windows-msvc/src/sirius/app/Release/sirius.exe'} else {'bin/windows-msvc/tests/backend/Release/sirius_backend_tests.exe'}
$loaderKey=$exeKey.Substring(0,$exeKey.LastIndexOf('/')+1)+'vulkan-1.dll';$loaderBinding=$manifest.files.PSObject.Properties[$loaderKey].Value
$loaderPath=Join-Path ([IO.Path]::GetDirectoryName($owner.test_executable)) 'vulkan-1.dll'
if (!$loaderBinding -or $owner.expected_staged_loader.sha256 -cne $loaderBinding.sha256 -or $owner.expected_staged_loader.bytes -ne $loaderBinding.bytes -or ![String]::Equals([IO.Path]::GetFullPath($owner.expected_staged_loader.path),[IO.Path]::GetFullPath($loaderPath),[StringComparison]::OrdinalIgnoreCase)) { throw 'Staged SDK loader authority differs' }
$births=New-Object 'System.Collections.Generic.List[object]'
foreach ($name in @('owner','test')) {
 $pidValue=if ($name -eq 'owner') {[int]$owner.owner_pid} else {[int]$owner.test_pid}
 $startValue=if ($name -eq 'owner') {$owner.owner_start_time_utc} else {$owner.test_start_time_utc}
 $p=Get-Process -Id $pidValue -ErrorAction SilentlyContinue
 $sameBirth=$false
 if ($p) { $sameBirth=($p.StartTime.ToUniversalTime().ToString('o') -eq $startValue) }
 $births.Add([ordered]@{kind=$name;pid=$pidValue;recorded_start_utc=$startValue;pid_currently_live=($null -ne $p);owned_birth_absent=(!$sameBirth)})
 if ($sameBirth) { throw 'Owned execution birth still live' }
}
$seals=New-Object 'System.Collections.Generic.List[object]'
$modulePath=Join-Path $Workspace ($RunName+'\driver-modules.json')
$actualModulesObserved=Test-Path -LiteralPath $modulePath
if ($owner.invocation_kind -eq 'Control' -and !$actualModulesObserved) { throw 'Original smoke requires actual loader and AMD provider capture' }
if ($actualModulesObserved) {
 $modules=Get-Content -Raw -LiteralPath $modulePath | ConvertFrom-Json
 foreach ($m in $modules) {
  $actual=(Get-FileHash -Algorithm SHA256 -LiteralPath $m.path).Hash.ToLowerInvariant()
  $actualBytes=(Get-Item -LiteralPath $m.path).Length
  $same=($actual -eq $m.sha256 -and $actualBytes -eq $m.bytes)
  $seals.Add([ordered]@{path=$m.path;expected_sha256=$m.sha256;post_sha256=$actual;expected_bytes=$m.bytes;post_bytes=$actualBytes;unchanged=$same})
  if (!$same) { throw 'Loaded provider module post-seal differs' }
  if ($m.name -ieq 'vulkan-1.dll' -and (![String]::Equals([IO.Path]::GetFullPath($m.path),[IO.Path]::GetFullPath($loaderPath),[StringComparison]::OrdinalIgnoreCase) -or $m.sha256 -cne $loaderBinding.sha256 -or $m.bytes -ne $loaderBinding.bytes)) { throw 'Observed loader is not the exact staged SDK loader' }
 }
 if ($seals.Count -ne 2) { throw 'Observed module receipt must contain both loader and AMD provider' }
}
$loader=Join-Path $env:WINDIR 'System32\vulkan-1.dll'
$installedLoader=[ordered]@{path=$loader;bytes=(Get-Item -LiteralPath $loader).Length;sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $loader).Hash.ToLowerInvariant();version=(Get-Item -LiteralPath $loader).VersionInfo.FileVersion;scope='Installed file input only; does not establish it was mapped by the inventory process.'}
$exeHash=(Get-FileHash -Algorithm SHA256 -LiteralPath $owner.test_executable).Hash.ToLowerInvariant()
if ($exeHash -ne $owner.test_executable_sha256) { throw 'Owned test executable post-seal differs' }
$processes=@(Get-CimInstance Win32_Process -OperationTimeoutSec 20)
$matches=New-Object 'System.Collections.Generic.List[object]';$compilers=New-Object 'System.Collections.Generic.List[object]';$hidden=0
foreach ($p in $processes) {
 if ([int]$p.ProcessId -eq $PID) { continue }
 if (!$p.CommandLine -or !$p.ExecutablePath) { $hidden++ }
 if ($p.Name -ieq 'csc.exe') { $compilers.Add([ordered]@{pid=$p.ProcessId;parent_pid=$p.ParentProcessId;creation_date=$p.CreationDate.ToUniversalTime().ToString('o');executable=$p.ExecutablePath;command_line=$p.CommandLine}) }
 if (($p.CommandLine -and $p.CommandLine.Contains($Workspace)) -or ($p.ExecutablePath -and $p.ExecutablePath.Contains($Workspace)) -or ($p.CommandLine -and $p.CommandLine.Contains($BuildRoot)) -or ($p.ExecutablePath -and $p.ExecutablePath.Contains($BuildRoot))) {
  $matches.Add([ordered]@{pid=$p.ProcessId;creation_date=$p.CreationDate.ToUniversalTime().ToString('o');executable=$p.ExecutablePath;command_line=$p.CommandLine})
 }
}
if ($compilers.Count) { throw 'Accessible CSharp compiler remains live' }
if ($matches.Count) { throw 'Accessible native consumer still references scratch workspace' }
$v=[ordered]@{schema_version=1;kind='terminal ownership/provider seal audit';run=$RunName;births=$births.ToArray();loaded_module_post_seals=$seals.ToArray();owned_executable_post_seal=[ordered]@{path=$owner.test_executable;sha256=$exeHash;unchanged=$true};owned_births_absent=$true;all_observed_module_seals_match=$true;actual_loaded_modules_observed=$actualModulesObserved;source_product_manifest_sha256=$manifestHash;native_build_or_mandatory_gate_admitted=$false;installed_loader_input=$installedLoader;consumer_scan=[ordered]@{processes=$processes.Count;hidden_field_records=$hidden;matches=$matches.ToArray();accessible_compilers=$compilers.ToArray()};scope='Diagnostic inventory or original smoke only; absent short-lived inventory module capture is explicitly incomplete loaded-provider evidence. Installed loader bytes do not replace it. Normal completion required. Accessible command/executable scan only; no handle enumeration, unconditional compiler-descendant guard ownership or all-consumer absence claim.'}
$v | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Join-Path $Workspace $Receipt) -Encoding UTF8
[pscustomobject]@{run=$RunName;owned_births_absent=$true;module_seals=$seals.Count;accessible_processes=$processes.Count;hidden_field_records=$hidden;consumer_matches=$matches.Count;accessible_compilers=$compilers.Count} | ConvertTo-Json -Compress
