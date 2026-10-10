param([Parameter(Mandatory=$true)][string]$Trial,[Parameter(Mandatory=$true)][string]$Receipt)
$ErrorActionPreference='Stop'
$known=New-Object 'System.Collections.Generic.List[object]'
foreach ($file in Get-ChildItem -LiteralPath (Join-Path $Trial 'native') -Filter owner.json -Recurse) {
 $owner=Get-Content -Raw -LiteralPath $file.FullName | ConvertFrom-Json
 foreach ($kind in @('owner','test')) {
  $number=if ($kind -eq 'owner') {[int]$owner.owner_pid} else {[int]$owner.test_pid}
  $birth=if ($kind -eq 'owner') {$owner.owner_start_time_utc} else {$owner.test_start_time_utc}
  $process=Get-Process -Id $number -ErrorAction SilentlyContinue
  $same=$false
  if ($process) { $same=($process.StartTime.ToUniversalTime().ToString('o') -eq $birth) }
  $known.Add([ordered]@{kind=$kind;pid=$number;recorded_start_utc=$birth;owned_birth_absent=(!$same)})
  if ($same) { throw 'Recorded native birth remains live' }
 }
}
if ($known.Count -ne 12) { throw 'Expected six completed native owner/test pairs' }
$targets=@('submission-fence-review','submission-fence-restoration','native-312bda3','native-c6b271d')
$processes=@(Get-CimInstance Win32_Process -OperationTimeoutSec 20)
$matches=New-Object 'System.Collections.Generic.List[object]'
$hidden=0
foreach ($process in $processes) {
 if ([int]$process.ProcessId -eq $PID) { continue }
 if (!$process.CommandLine -or !$process.ExecutablePath) { $hidden++ }
 foreach ($target in $targets) {
  if (($process.CommandLine -and $process.CommandLine.IndexOf($target,[StringComparison]::OrdinalIgnoreCase) -ge 0) -or ($process.ExecutablePath -and $process.ExecutablePath.IndexOf($target,[StringComparison]::OrdinalIgnoreCase) -ge 0)) {
   $matches.Add([ordered]@{pid=$process.ProcessId;target=$target;executable=$process.ExecutablePath;command_line=$process.CommandLine})
  }
 }
}
if ($matches.Count) { throw 'Accessible scoped native consumer remains live' }
$value=[ordered]@{checked_utc=[DateTime]::UtcNow.ToString('o');pass=$true;known_births=$known.ToArray();targets=$targets;process_count=$processes.Count;hidden_field_records=$hidden;matches=$matches.ToArray();scope='Fresh recorded-birth and accessible command/executable scan only. Self observer excluded; inaccessible fields and handles/global/scheduled consumers unverified.'}
$value | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $Receipt -Encoding UTF8
[pscustomobject]@{pass=$true;known_births=$known.Count;process_count=$processes.Count;hidden_field_records=$hidden;matches=$matches.Count} | ConvertTo-Json -Compress
