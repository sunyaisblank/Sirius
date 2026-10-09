param([int]$TestId,[string]$StartTime)
$p=Get-Process -Id $TestId -ErrorAction SilentlyContinue
$absent=(!$p) -or ($p.StartTime.ToUniversalTime().ToString('o') -ne $StartTime)
[ordered]@{test_pid=$TestId;original_start_utc=$StartTime;owned_process_absent=$absent;checked_utc=[DateTime]::UtcNow.ToString('o');scope='Known PID and creation time only; no process scan or termination'} | ConvertTo-Json
if (!$absent) { exit 1 }
