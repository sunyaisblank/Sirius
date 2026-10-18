param([Parameter(Mandatory=$true)][string]$Bootstrap)
$ErrorActionPreference='Stop'
[Console]::OutputEncoding=New-Object System.Text.UTF8Encoding($false)
$b=Get-Content -Raw -LiteralPath $Bootstrap | ConvertFrom-Json
if ($b.schema -ne 'sirius-july-nt-snapshot-name-bootstrap-v1' -or [int]$b.owner_pid -le 0) { throw 'Wrong owned bootstrap' }
$saved=([DateTimeOffset]::Parse($b.owner_start_utc)).UtcDateTime.Ticks
$r=[ordered]@{schema='sirius-july-nt-snapshot-name-terminal-v1';owner_pid=[int]$b.owner_pid;saved_start_utc=$b.owner_start_utc;matched_live_birth=$false;stopped_owned_birth=$false;owned_birth_absent=$false;scope='Exact task PowerShell birth only; no global process census or VSS absence claim.'}
$p=$null
try {
    try { $p=[Diagnostics.Process]::GetProcessById([int]$b.owner_pid);$held=$p.Handle } catch [ArgumentException] { $r.owned_birth_absent=$true }
    if ($p) {
        $actual=$p.StartTime.ToUniversalTime().Ticks
        if ($actual -ne $saved) { $r.owned_birth_absent=$true } else {
            $r.matched_live_birth=$true
            $p.Kill()
            if (-not $p.WaitForExit(8000)) { throw 'Owned PowerShell birth did not retire' }
            $r.stopped_owned_birth=$true;$r.owned_birth_absent=$true
        }
    }
} finally { if ($p) { $p.Dispose() } }
$r | ConvertTo-Json -Compress
if (-not $r.owned_birth_absent) { exit 2 }
