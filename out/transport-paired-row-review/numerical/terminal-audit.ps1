param([Parameter(Mandatory=$true)][string]$Workspace,[string]$RunName='run',[string]$Receipt='terminal-census-audit.json')
$ErrorActionPreference='Stop'
$owner=Get-Content -Raw -LiteralPath (Join-Path $Workspace ($RunName+'\owner.json')) | ConvertFrom-Json
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
$modules=Get-Content -Raw -LiteralPath (Join-Path $Workspace ($RunName+'\driver-modules.json')) | ConvertFrom-Json
foreach ($m in $modules) {
 $actual=(Get-FileHash -Algorithm SHA256 -LiteralPath $m.path).Hash.ToLowerInvariant()
 $actualBytes=(Get-Item -LiteralPath $m.path).Length
 $same=($actual -eq $m.sha256 -and $actualBytes -eq $m.bytes)
 $seals.Add([ordered]@{path=$m.path;expected_sha256=$m.sha256;post_sha256=$actual;expected_bytes=$m.bytes;post_bytes=$actualBytes;unchanged=$same})
 if (!$same) { throw 'Loaded provider module post-seal differs' }
}
if ($seals.Count -ne 2) { throw 'Expected exactly two actual loaded module seals' }
$exeHash=(Get-FileHash -Algorithm SHA256 -LiteralPath $owner.test_executable).Hash.ToLowerInvariant()
if ($exeHash -ne $owner.test_executable_sha256) { throw 'Owned PowerShell executable post-seal differs' }
$processes=@(Get-CimInstance Win32_Process)
$matches=New-Object 'System.Collections.Generic.List[object]';$compilers=New-Object 'System.Collections.Generic.List[object]';$hidden=0
foreach ($p in $processes) {
 if ([int]$p.ProcessId -eq $PID) { continue }
 if (!$p.CommandLine -or !$p.ExecutablePath) { $hidden++ }
 if ($p.Name -ieq 'csc.exe') { $compilers.Add([ordered]@{pid=$p.ProcessId;parent_pid=$p.ParentProcessId;creation_date=$p.CreationDate.ToUniversalTime().ToString('o');executable=$p.ExecutablePath;command_line=$p.CommandLine}) }
 if (($p.CommandLine -and $p.CommandLine.Contains($Workspace)) -or ($p.ExecutablePath -and $p.ExecutablePath.Contains($Workspace))) {
  $matches.Add([ordered]@{pid=$p.ProcessId;creation_date=$p.CreationDate.ToUniversalTime().ToString('o');executable=$p.ExecutablePath;command_line=$p.CommandLine})
 }
}
if ($matches.Count) { throw 'Accessible native consumer still references scratch workspace' }
$v=[ordered]@{schema_version=1;kind='terminal ownership/provider seal audit';run=$RunName;births=$births.ToArray();loaded_module_post_seals=$seals.ToArray();owned_executable_post_seal=[ordered]@{path=$owner.test_executable;sha256=$exeHash;unchanged=$true};owned_births_absent=$true;all_two_module_seals_match=$true;consumer_scan=[ordered]@{processes=$processes.Count;hidden_field_records=$hidden;matches=$matches.ToArray();accessible_compilers=$compilers.ToArray()};scope='Normal completion required. Accessible command/executable scan only; no handle enumeration, unconditional compiler-descendant guard ownership or all-consumer absence claim.'}
$v | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Join-Path $Workspace $Receipt) -Encoding UTF8
[pscustomobject]@{run=$RunName;owned_births_absent=$true;module_seals=$seals.Count;accessible_processes=$processes.Count;hidden_field_records=$hidden;consumer_matches=$matches.Count;accessible_compilers=$compilers.Count} | ConvertTo-Json -Compress
