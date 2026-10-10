param([Parameter(Mandatory=$true)][string]$Workspace,[Parameter(Mandatory=$true)][string]$BaselineStage,[Parameter(Mandatory=$true)][string]$CandidateStage,[Parameter(Mandatory=$true)][string]$Receipt)
$ErrorActionPreference='Stop'
$keys=New-Object 'System.Collections.Generic.HashSet[string]'
$births=New-Object 'System.Collections.Generic.List[object]'
function Check-Birth([int]$ProcessId,[string]$Start,[string]$Kind) {
 if (!$keys.Add(([string]$ProcessId+'|'+$Start))) { return }
 $p=Get-Process -Id $ProcessId -ErrorAction SilentlyContinue
 $same=$false
 if ($p) { $same=($p.StartTime.ToUniversalTime().ToString('o') -eq $Start) }
 if ($same) { throw 'Recorded native execution or auditor birth remains live' }
 $births.Add([ordered]@{kind=$Kind;pid=$ProcessId;recorded_start_utc=$Start;owned_birth_absent=$true})
}
$owners=@(Get-ChildItem -LiteralPath (Join-Path $Workspace 'native') -Filter 'owner.json' -Recurse -File)
if ($owners.Count -ne 15) { throw 'Expected fifteen completed native owners' }
$exes=@{};$modules=@{}
foreach ($path in $owners) {
 $o=Get-Content -Raw -LiteralPath $path.FullName | ConvertFrom-Json
 if ($o.status -ne 'completed' -or $o.returncode -ne 0 -or $o.outer_stop -or $o.cleanup_errors.Count) { throw 'Native owner did not complete normally' }
 Check-Birth $o.owner_pid $o.owner_start_time_utc 'owner'
 Check-Birth $o.test_pid $o.test_start_time_utc 'test'
 $exes[$o.test_executable]=$o.test_executable_sha256
 foreach ($m in (Get-Content -Raw -LiteralPath (Join-Path $path.DirectoryName 'driver-modules.json') | ConvertFrom-Json)) { $modules[$m.path]=$m }
}
$bootstraps=@(Get-ChildItem -LiteralPath (Join-Path $Workspace 'temp') -Filter 'bootstrap.json' -Recurse -File)
if ($bootstraps.Count -ne 30) { throw 'Expected thirty historical bootstrap receipts' }
foreach ($path in $bootstraps) {
 $b=Get-Content -Raw -LiteralPath $path.FullName | ConvertFrom-Json
 Check-Birth $b.owner_pid $b.owner_start_time_utc 'bootstrap-owner-or-auditor'
}
if ($births.Count -ne 45 -or $exes.Count -ne 2 -or $modules.Count -ne 2) { throw 'Native birth/product/provider closure differs' }
$seals=New-Object 'System.Collections.Generic.List[object]'
foreach ($path in $exes.Keys) {
 $actual=(Get-FileHash -Algorithm SHA256 -LiteralPath $path).Hash.ToLowerInvariant()
 if ($actual -ne $exes[$path]) { throw 'Native executable post-seal differs' }
 $seals.Add([ordered]@{path=$path;sha256=$actual;kind='native-test-executable';unchanged=$true})
}
foreach ($path in $modules.Keys) {
 $actual=(Get-FileHash -Algorithm SHA256 -LiteralPath $path).Hash.ToLowerInvariant();$m=$modules[$path]
 if ($actual -ne $m.sha256 -or (Get-Item -LiteralPath $path).Length -ne $m.bytes) { throw 'Native provider post-seal differs' }
 $seals.Add([ordered]@{path=$path;sha256=$actual;bytes=$m.bytes;kind='loaded-provider';unchanged=$true})
}
$processes=@(Get-CimInstance Win32_Process -OperationTimeoutSec 20);$matches=@();$compilers=@();$hidden=0
foreach ($p in $processes) {
 if ([int]$p.ProcessId -eq $PID) { continue }
 if (!$p.CommandLine -or !$p.ExecutablePath) { $hidden++ }
 if ($p.Name -ieq 'csc.exe') { $compilers+=@{pid=$p.ProcessId;creation_date=$p.CreationDate.ToUniversalTime().ToString('o');executable=$p.ExecutablePath;command_line=$p.CommandLine} }
 foreach ($root in @($Workspace,$BaselineStage,$CandidateStage)) {
  if (($p.CommandLine -and $p.CommandLine.IndexOf($root,[StringComparison]::OrdinalIgnoreCase) -ge 0) -or ($p.ExecutablePath -and $p.ExecutablePath.IndexOf($root,[StringComparison]::OrdinalIgnoreCase) -ge 0)) { $matches+=@{pid=$p.ProcessId;creation_date=$p.CreationDate.ToUniversalTime().ToString('o');executable=$p.ExecutablePath;command_line=$p.CommandLine};break }
 }
}
if ($matches.Count -or $compilers.Count) { throw 'Accessible task consumer or CSharp compiler remains live' }
@{pass=$true;births=$births.ToArray();post_seals=$seals.ToArray();census=@{processes=$processes.Count;hidden_field_records=$hidden;matches=$matches;accessible_compilers=$compilers};scope='Fresh read-only exact recorded Windows birth absence and accessible command/executable reference census for finished trial and both native stages. No handle enumeration or hidden/global/unconditional descendant absence claim.'} | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath $Receipt -Encoding UTF8
[pscustomobject]@{pass=$true;recorded_births_absent=$births.Count;provider_and_product_seals=$seals.Count;processes=$processes.Count;hidden_field_records=$hidden}|ConvertTo-Json -Compress
