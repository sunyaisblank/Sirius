param([Parameter(Mandatory=$true)][string]$Workspace)
$ErrorActionPreference='Stop'
$self=Get-Process -Id $PID
$handle=$self.Handle
@{owner_pid=$PID;owner_start_time_utc=$self.StartTime.ToUniversalTime().ToString('o');expected_test_executable='';scope='Read-only obsolete-stage consumer observer; no compiler or science launch'} | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $Workspace 'temp\bootstrap.json') -Encoding UTF8
$roots=Get-Content -Raw -LiteralPath (Join-Path $Workspace 'windows-roots.json') | ConvertFrom-Json
$processes=@(Get-CimInstance Win32_Process -OperationTimeoutSec 20)
$matches=New-Object 'System.Collections.Generic.List[object]'
$hidden=0
foreach ($p in $processes) {
 if ([int]$p.ProcessId -eq $PID) { continue }
 if (!$p.CommandLine -or !$p.ExecutablePath) { $hidden++ }
 foreach ($root in $roots) {
  if (($p.CommandLine -and $p.CommandLine.IndexOf($root,[StringComparison]::OrdinalIgnoreCase) -ge 0) -or ($p.ExecutablePath -and $p.ExecutablePath.IndexOf($root,[StringComparison]::OrdinalIgnoreCase) -ge 0)) {
   $matches.Add([ordered]@{root=$root;pid=$p.ProcessId;executable=$p.ExecutablePath;command_line=$p.CommandLine})
  }
 }
}
@{pass=($matches.Count -eq 0);roots=[string[]]$roots;processes=$processes.Count;hidden_field_records=$hidden;matches=$matches.ToArray();scope='Accessible command/executable census only; no hidden-consumer, handle or global absence claim'} | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Workspace 'windows-consumers.json') -Encoding UTF8
if ($matches.Count) { throw 'Accessible consumer of obsolete stage remains live' }
