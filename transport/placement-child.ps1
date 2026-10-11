param([Parameter(Mandatory=$true)][string]$Stage,[Parameter(Mandatory=$true)][string]$OwnerTemporaryDirectory)
$ErrorActionPreference='Stop'
[Console]::OutputEncoding=[Text.UTF8Encoding]::new($false)
$self=Get-Process -Id $PID;$held=$self.Handle;$birth=$self.StartTime.ToUniversalTime().ToString('o')
try {
 $nativeTemp=New-Item -ItemType Directory -Path (Join-Path $Stage 'native-temp')
 $env:TEMP=$nativeTemp.FullName;$env:TMP=$env:TEMP
 $bootstrap=[ordered]@{child_pid=$PID;child_start_time_utc=$birth;scope='Exact child birth held before compiler/provider execution'}
 [IO.File]::WriteAllText((Join-Path $Stage 'child-bootstrap.json'),($bootstrap|ConvertTo-Json),[Text.UTF8Encoding]::new($false))
 $owner=Get-Content -Raw -LiteralPath (Join-Path $OwnerTemporaryDirectory 'bootstrap.json')|ConvertFrom-Json
 $clock=[Diagnostics.Stopwatch]::StartNew()
 while(!(Test-Path -LiteralPath (Join-Path $Stage 'START.json'))) {
  if($clock.Elapsed.TotalSeconds -gt 20){throw 'Parent START barrier exceeded20s'}
  $parent=Get-Process -Id $owner.owner_pid -ErrorAction SilentlyContinue
  if(!$parent -or $parent.StartTime.ToUniversalTime().ToString('o') -ne $owner.owner_start_time_utc){throw 'Recorded parent absent before START'}
  Start-Sleep -Milliseconds 20
 }
 $start=Get-Content -Raw -LiteralPath (Join-Path $Stage 'START.json')|ConvertFrom-Json
 if($start.child_pid -ne $PID -or $start.child_start_time_utc -ne $birth -or $start.owner_pid -ne $owner.owner_pid -or $start.owner_start_time_utc -ne $owner.owner_start_time_utc){throw 'START does not bind exact owner/child births'}
 Add-Type -TypeDefinition ([IO.File]::ReadAllText((Join-Path $Stage 'native_clock_placement.cs')))
 $payloads=(New-Item -ItemType Directory -Path (Join-Path $Stage 'query-payloads')).FullName
 $result=[NativeClockPlacement]::Observe('C:\Windows\System32\vulkan-1.dll',(Join-Path $Stage 'retained_transport.spv'),(Join-Path $Stage 'retained_endpoint.spv'),$payloads,[int]$owner.owner_pid,[string]$owner.owner_start_time_utc)
 [IO.File]::WriteAllText((Join-Path $Stage 'clock-placement.json'),($result|ConvertTo-Json -Depth 10),[Text.UTF8Encoding]::new($false))
 [IO.File]::WriteAllText((Join-Path $Stage 'run/driver-modules.json'),(ConvertTo-Json -InputObject @($result['loaded_modules']) -Depth 8),[Text.UTF8Encoding]::new($false))
 [Console]::WriteLine('Native clock placement query completed')
 [Environment]::Exit(0)
}catch{[Console]::Error.WriteLine($_.Exception.ToString());[Environment]::Exit(1)}
