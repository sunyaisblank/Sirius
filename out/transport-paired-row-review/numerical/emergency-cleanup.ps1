param([Parameter(Mandatory=$true)][string]$TemporaryDirectory,[Parameter(Mandatory=$true)][string]$Output,[Parameter(Mandatory=$true)][string]$Receipt,[Parameter(Mandatory=$true)][string]$TaskWorkspace)
$ErrorActionPreference='Stop'
$record=[ordered]@{kind='invalidated-native-owner-emergency-cleanup';helper_pid=$PID;helper_start_time_utc=(Get-Process -Id $PID).StartTime.ToUniversalTime().ToString('o');targets=@();errors=@();scope='Exact recorded owner/test births and accessible task-workspace compiler command lines only; no hidden-handle or unconditional compiler-descendant containment claim.'}
$errors=New-Object 'System.Collections.Generic.List[string]'
$targets=New-Object 'System.Collections.Generic.List[object]'
try {
 $bootstrapPath=Join-Path $TemporaryDirectory 'bootstrap.json'
 if (!(Test-Path -LiteralPath $bootstrapPath)) { throw 'Early native owner birth unavailable; do not resume dependent device work' }
 $bootstrap=Get-Content -Raw -LiteralPath $bootstrapPath | ConvertFrom-Json
 $ownerPath=Join-Path $Output 'owner.json'
 try {
 if (Test-Path -LiteralPath $ownerPath) {
  $owner=Get-Content -Raw -LiteralPath $ownerPath | ConvertFrom-Json
  if ($owner.test_pid -and $owner.test_start_time_utc) { $targets.Add([ordered]@{kind='test';pid=[int]$owner.test_pid;birth=$owner.test_start_time_utc}) }
 }
 } catch { $errors.Add('optional-owner-receipt: '+$_.Exception.Message) }
 try {
 $rootProcess=Get-Process -Id ([int]$bootstrap.owner_pid) -ErrorAction SilentlyContinue
 if ($rootProcess) { $rootHandle=$rootProcess.Handle }
 if ($rootProcess -and $rootProcess.StartTime.ToUniversalTime().ToString('o') -ne $bootstrap.owner_start_time_utc) { throw 'Native owner PID was reused; refuse descendant targeting' }
 $initialProcesses=@(Get-CimInstance Win32_Process -OperationTimeoutSec 20)
 foreach ($item in $initialProcesses) {
  if ([int]$item.ParentProcessId -ne [int]$bootstrap.owner_pid -or $item.CreationDate.ToUniversalTime() -lt [DateTime]$bootstrap.owner_start_time_utc) { continue }
  $isTest=($item.ExecutablePath -and [String]::Equals($item.ExecutablePath,$bootstrap.expected_test_executable,[StringComparison]::OrdinalIgnoreCase))
  $isCompiler=($item.Name -ieq 'csc.exe' -and $item.CommandLine -and $item.CommandLine.IndexOf($TaskWorkspace,[StringComparison]::OrdinalIgnoreCase) -ge 0)
  if ($isTest -or $isCompiler) {
   if (@($targets | Where-Object {$_.pid -eq [int]$item.ProcessId}).Count -eq 0) {
    $discovered=Get-Process -Id $item.ProcessId -ErrorAction SilentlyContinue
    if ($discovered) {
     $discoveredHandle=$discovered.Handle
     $actualBirth=$discovered.StartTime.ToUniversalTime()
     if ([Math]::Truncate([decimal]$actualBirth.Ticks / 10) -ne [Math]::Truncate([decimal]$item.CreationDate.ToUniversalTime().Ticks / 10)) { throw 'Discovered PID birth differs at CIM microsecond precision' }
     $targets.Add([ordered]@{kind=if ($isTest) {'discovered-test'} else {'accessible-task-compiler'};pid=[int]$item.ProcessId;birth=$actualBirth.ToString('o');cim_birth=$item.CreationDate.ToUniversalTime().ToString('o')})
    }
   }
  } else { throw 'Owned parent has an unidentified accessible child; do not resume dependent work' }
 }
 } catch { $errors.Add('descendant-discovery: '+$_.Exception.Message) }
 $targets.Add([ordered]@{kind='owner';pid=[int]$bootstrap.owner_pid;birth=$bootstrap.owner_start_time_utc})
 foreach ($target in $targets) {
  try {
  $p=Get-Process -Id $target.pid -ErrorAction SilentlyContinue
  $target.same_birth_live=$false
  if ($p) { $heldHandle=$p.Handle }
  if ($p -and $p.StartTime.ToUniversalTime().ToString('o') -eq $target.birth) {
   $target.same_birth_live=$true
   $p.Kill()
   if (!$p.WaitForExit(5000)) { throw 'Recorded owned birth did not terminate' }
  }
  $p=Get-Process -Id $target.pid -ErrorAction SilentlyContinue
  $target.owned_birth_absent=(!$p -or $p.StartTime.ToUniversalTime().ToString('o') -ne $target.birth)
  if (!$target.owned_birth_absent) { throw 'Recorded owned birth remains live' }
  } catch { $errors.Add('recorded-birth-cleanup: '+$_.Exception.Message) }
 }
 $processes=@(Get-CimInstance Win32_Process -OperationTimeoutSec 20)
 $hidden=0;$compilerStops=New-Object 'System.Collections.Generic.List[object]'
 foreach ($item in $processes) {
  if (!$item.CommandLine -or !$item.ExecutablePath) { $hidden++ }
  if ($item.Name -ieq 'csc.exe' -and $item.CommandLine -and $item.CommandLine.IndexOf($TaskWorkspace,[StringComparison]::OrdinalIgnoreCase) -ge 0) {
   $p=Get-Process -Id $item.ProcessId -ErrorAction SilentlyContinue
   if ($p) {
    $heldHandle=$p.Handle
    $actualBirth=$p.StartTime.ToUniversalTime()
    if ([Math]::Truncate([decimal]$actualBirth.Ticks / 10) -ne [Math]::Truncate([decimal]$item.CreationDate.ToUniversalTime().Ticks / 10)) { throw 'Task compiler birth differs at CIM microsecond precision' }
    $p.Kill();if (!$p.WaitForExit(5000)) { throw 'Accessible task compiler did not terminate' }
    $still=Get-Process -Id $item.ProcessId -ErrorAction SilentlyContinue
    if ($still -and $still.StartTime.ToUniversalTime().ToString('o') -eq $actualBirth.ToString('o')) { throw 'Task compiler owned birth remains live' }
    $compilerStops.Add([ordered]@{pid=$item.ProcessId;birth=$actualBirth.ToString('o');cim_birth=$item.CreationDate.ToUniversalTime().ToString('o');command=$item.CommandLine;owned_birth_absent=$true})
   }
  }
 }
 $record.visible_processes=$processes.Count;$record.hidden_field_records=$hidden;$record.accessible_task_compiler_stops=$compilerStops.ToArray()
} catch { $errors.Add($_.Exception.Message) }
$record.targets=$targets.ToArray();$record.errors=$errors.ToArray();$record.recorded_births_absent=($targets.Count -gt 0 -and @($targets | Where-Object {!$_.owned_birth_absent}).Count -eq 0)
$record | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath $Receipt -Encoding UTF8
if ($errors.Count -or !$record.recorded_births_absent) { exit 1 }
