param(
    [Parameter(Mandatory=$true)][string]$Inputs,
    [Parameter(Mandatory=$true)][string]$Receipt
)
$ErrorActionPreference='Stop'
$data=Get-Content -Raw -LiteralPath $Inputs | ConvertFrom-Json
$known=@()
foreach($birth in @($data.windows_births)){
    $process=$null
    try { $process=Get-Process -Id ([int]$birth.pid) -ErrorAction Stop }
    catch { if($_.FullyQualifiedErrorId -notlike 'NoProcessFoundForGivenId*'){throw} }
    $actual=$null
    if($null -ne $process){
        $handle=$process.Handle
        $actual=$process.StartTime.ToUniversalTime().ToString('o')
        if($process.StartTime.ToUniversalTime().Ticks -eq ([DateTime]::Parse([string]$birth.start_utc)).ToUniversalTime().Ticks){
            throw ('Recorded native birth remains: '+$birth.pid)
        }
    }
    $known+=@{pid=[int]$birth.pid;saved_start_utc=[string]$birth.start_utc;current_start_utc=$actual;recorded_birth_absent=$true}
}
$pins=@()
foreach($pin in @($data.provider_pins)){
    $file=Get-Item -LiteralPath ([string]$pin.windows_path)
    $hash=(Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    if($file.Length -ne [long]$pin.bytes -or $hash -ne [string]$pin.sha256){throw ('Provider post seal mismatch: '+$pin.windows_path)}
    $pins+=@{path=$file.FullName;bytes=$file.Length;sha256=$hash;matches=$true}
}
$matches=@();$compilers=@();$hidden=0;$count=0
foreach($process in @(Get-CimInstance Win32_Process -OperationTimeoutSec 20)){
    $count++
    if($process.ProcessId -eq $PID){continue}
    $command=[string]$process.CommandLine;$executable=[string]$process.ExecutablePath
    if([string]::IsNullOrEmpty($command) -or [string]::IsNullOrEmpty($executable)){$hidden++}
    $fields=@()
    foreach($target in @($data.windows_targets)){
        if($command.IndexOf([string]$target,[StringComparison]::OrdinalIgnoreCase) -ge 0){$fields+='command'}
        if($executable.IndexOf([string]$target,[StringComparison]::OrdinalIgnoreCase) -ge 0){$fields+='executable'}
    }
    if($fields.Count -gt 0){$matches+=@{pid=$process.ProcessId;name=$process.Name;fields=$fields}}
    if([string]$process.Name -ieq 'csc.exe'){$compilers+=@{pid=$process.ProcessId;command=$command;executable=$executable}}
}
if($matches.Count -ne 0 -or $compilers.Count -ne 0){throw 'Accessible removal-target consumer/compiler remains'}
@{pass=$true;utc=[DateTime]::UtcNow.ToString('o');known_births=$known;known_births_verified=$known.Count;provider_pins=$pins;processes=$count;hidden_field_records=$hidden;matches=$matches;accessible_compilers=$compilers;scope='Exact recorded Windows births with acquired handles before live birth checks, three provider post seals and bounded accessible command/executable census. No scientific execution/Add-Type, hidden-field/global-handle or unconditional compiler containment claim.'}|ConvertTo-Json -Depth 10|Set-Content -LiteralPath $Receipt -Encoding UTF8
