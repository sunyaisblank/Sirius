$ErrorActionPreference = 'Stop'
foreach ($item in @(@{path='\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\retained-qualification\native-backend-125a890\observation\process.json'; name='sirius_backend_tests'}, @{path='\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\retained-qualification\native-backend-125a890\observation\inventory-process.json'; name='sirius'}, @{path='\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\retained-qualification\native-backend-125a890\observation\supervisor-process.json'; name='powershell'})) {
    if (-not (Test-Path -LiteralPath $item.path)) { continue }
    $state = Get-Content -Raw -LiteralPath $item.path | ConvertFrom-Json
    $p = Get-Process -Id $state.pid -ErrorAction SilentlyContinue
    if ($null -eq $p) { continue }
    $null = $p.Handle
    if ($p.StartTime.ToUniversalTime().Ticks -ne [int64]$state.start_time_ticks) { throw 'PID recycled; refusal to terminate unrelated process' }
    if ($p.ProcessName -ne $item.name -or $p.MainModule.FileName -ne $state.executable) { throw 'Exact native process path mismatch' }
    $p.Kill()
    if (-not $p.WaitForExit(5000)) { throw 'Exact native process did not terminate' }
}
exit 0
