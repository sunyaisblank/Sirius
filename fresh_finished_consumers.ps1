param([Parameter(Mandatory=$true)][string]$Work,[Parameter(Mandatory=$true)][string]$Births,[Parameter(Mandatory=$true)][string]$Output)
$ErrorActionPreference='Stop'
if (Test-Path -LiteralPath $Output) {throw 'Fresh census output exists'}
$Output=(New-Item -ItemType Directory -Path $Output).FullName
$temp=(New-Item -ItemType Directory -Path (Join-Path $Output 'temp')).FullName
$self=Get-Process -Id $PID;$held=$self.Handle
$bootstrap=[ordered]@{owner_pid=$PID;owner_start_time_utc=$self.StartTime.ToUniversalTime().ToString('o');expected_test_executable='';scope='Read-only census owner; no compiler, child tool or scientific execution'}
$bootstrap | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $temp 'bootstrap.json') -Encoding UTF8
$record=[ordered]@{owner_pid=$PID;owner_start_time_utc=$bootstrap.owner_start_time_utc;pass=$false;errors=@();known_births=@();matches=@();scope='Fresh recorded births and finite accessible Windows consumer fields only; hidden fields, handles and scheduled/global consumers remain unverified'}
try {
 $known=Get-Content -Raw -LiteralPath $Births | ConvertFrom-Json
 if (!$known.windows.Count) {throw 'No recorded Windows births'}
 foreach ($entry in $known.windows) {
  $process=Get-Process -Id ([int]$entry.pid) -ErrorAction SilentlyContinue
  $same=$false
  if ($process) {$handle=$process.Handle;$same=($process.StartTime.ToUniversalTime().Ticks -eq ([DateTime]$entry.birth).ToUniversalTime().Ticks)}
  $record.known_births+=@([ordered]@{pid=$entry.pid;birth=$entry.birth;owned_birth_absent=(!$same)})
  if ($same) {throw 'Recorded Windows birth remains live'}
 }
 $targets=@('native-spm-cost','native-b81e061-executed-cost','rdts-2026-10-02-1872')
 $processes=@(Get-CimInstance Win32_Process -OperationTimeoutSec 20)
 $record.targets=$targets;$record.process_count=$processes.Count;$record.hidden_field_records=0
 foreach ($process in $processes) {
  if ([int]$process.ProcessId -eq $PID) {continue}
  if (!$process.CommandLine -or !$process.ExecutablePath) {$record.hidden_field_records++}
  $match=($process.Name -match '^(sirius.*|RadeonDeveloper.*|RadeonGPUProfiler.*)\.exe$')
  foreach ($target in $targets) {
   if (($process.CommandLine -and $process.CommandLine.IndexOf($target,[StringComparison]::OrdinalIgnoreCase) -ge 0) -or ($process.ExecutablePath -and $process.ExecutablePath.IndexOf($target,[StringComparison]::OrdinalIgnoreCase) -ge 0)) {$match=$true}
  }
  if ($match) {$record.matches+=@([ordered]@{pid=$process.ProcessId;name=$process.Name;executable=$process.ExecutablePath;command_line=$process.CommandLine})}
 }
 if ($record.matches.Count) {throw 'Accessible Sirius/developer or scoped consumer remains live'}
 $record.pass=$true
} catch {$record.errors+=@($_.Exception.Message)}
$record.checked_utc=[DateTime]::UtcNow.ToString('o')
$record | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
[pscustomobject]@{pass=$record.pass;known_births=$record.known_births.Count;matches=$record.matches.Count;hidden_field_records=$record.hidden_field_records}|ConvertTo-Json -Compress
if (!$record.pass) {exit 1}
