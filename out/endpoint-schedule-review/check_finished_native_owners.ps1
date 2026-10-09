param([Parameter(Mandatory=$true)][string]$EvidenceRoot)
$ErrorActionPreference='Stop'
$records=@()
foreach ($file in Get-ChildItem -LiteralPath $EvidenceRoot -Filter owner.json -Recurse) {
 $owner=Get-Content -Raw -LiteralPath $file.FullName | ConvertFrom-Json
 if (!$owner.test_pid -or !$owner.test_start_time_utc) { continue }
 $p=Get-Process -Id $owner.test_pid -ErrorAction SilentlyContinue
 $absent=(!$p) -or ($p.StartTime.ToUniversalTime().ToString('o') -ne $owner.test_start_time_utc)
 $records+= [pscustomobject]@{control=$file.Directory.Name;test_pid=$owner.test_pid;recorded_start_utc=$owner.test_start_time_utc;owned_process_absent=$absent}
 if (!$absent) { throw 'Known recorded native child is still alive' }
}
[ordered]@{checked_utc=[DateTime]::UtcNow.ToString('o');scope='Known recorded child PIDs and exact creation times only; no scan or termination';records=$records} | ConvertTo-Json -Depth 6
