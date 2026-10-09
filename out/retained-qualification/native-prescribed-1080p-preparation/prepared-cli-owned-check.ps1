$ErrorActionPreference = 'Stop'
$recover = $false
$records = @()
try {
    foreach ($item in @(@{path='\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\retained-qualification\native-prescribed-1080p-preparation\cli-observation\process.json'; name='sirius'}, @{path='\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\retained-qualification\native-prescribed-1080p-preparation\cli-observation\inventory-process.json'; name='sirius'}, @{path='\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\retained-qualification\native-prescribed-1080p-preparation\cli-observation\supervisor-process.json'; name='powershell'})) {
        if (-not (Test-Path -LiteralPath $item.path)) { continue }
        $state = Get-Content -Raw -LiteralPath $item.path | ConvertFrom-Json
        $p = Get-Process -Id $state.pid -ErrorAction SilentlyContinue
        if ($null -eq $p) {
            $records += @{pid=$state.pid; start_time_ticks=$state.start_time_ticks; state='absent'}
            continue
        }
        $null = $p.Handle
        if ($p.StartTime.ToUniversalTime().Ticks -ne [int64]$state.start_time_ticks) { throw 'PID recycled; refusal to terminate unrelated process' }
        if ($p.ProcessName -ne $item.name -or $p.MainModule.FileName -ne $state.executable) { throw 'Exact native process path mismatch' }
        if (-not $recover) { throw 'An exact owned native process remains live' }
        $p.Kill()
        if (-not $p.WaitForExit(5000)) { throw 'Exact owned native process did not terminate in reap grace' }
        $records += @{pid=$state.pid; start_time_ticks=$state.start_time_ticks; state='owned-killed-and-reaped'}
        $p.Dispose()
    }
    [IO.File]::WriteAllText('\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\retained-qualification\native-prescribed-1080p-preparation\cli-observation\owned-process-check.json', (ConvertTo-Json -InputObject $records -Depth 4), [Text.UTF8Encoding]::new($false))
    exit 0
} catch {
    [Console]::Error.WriteLine($_.Exception.ToString())
    exit 1
}
