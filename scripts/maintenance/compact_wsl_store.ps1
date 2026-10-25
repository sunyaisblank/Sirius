# Run from Windows PowerShell. InspectOnly never stops a distribution.
# For an offline compaction, load this script into a scriptblock before
# stopping Ubuntu, because its checkout path is inside that distribution.
[CmdletBinding()]
param(
    [string] $DistributionName = 'Ubuntu',
    [switch] $InspectOnly,
    [switch] $StopDistribution
)
$ErrorActionPreference = 'Stop'
$entries = @(Get-ChildItem 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Lxss' |
    Get-ItemProperty | Where-Object { $_.DistributionName -eq $DistributionName })
if ($entries.Count -ne 1) { throw 'Expected one registered distribution.' }
$entry = $entries[0]
$filename = if ($entry.VhdFileName) { $entry.VhdFileName } else { 'ext4.vhdx' }
$path = Join-Path $entry.BasePath $filename
$file = Get-Item -LiteralPath $path -Force
if ($file.Extension -ne '.vhdx') { throw 'Expected the registered VHDX file.' }
function Get-RunningDistributions {
    @((& wsl.exe --list --running --quiet) | ForEach-Object { ($_ -replace "`0", '').Trim() } |
        Where-Object { $_ })
}
$running = $DistributionName -in (Get-RunningDistributions)
$drive = $file.PSDrive.Name
$beforeFree = (Get-PSDrive -Name $drive).Free
$report = [ordered]@{
    Distribution = $DistributionName
    Path = $file.FullName
    Running = $running
    Attributes = [string]$file.Attributes
    LengthBefore = $file.Length
    HostFreeBefore = $beforeFree
    OptimizeVhdAvailable = [bool](Get-Command Optimize-VHD -ErrorAction SilentlyContinue)
    WouldInterrupt = 'All work in the selected Ubuntu distribution, including Codex.'
    Completed = $false
}
if ($InspectOnly) { $report | ConvertTo-Json; return }
if (-not $report.OptimizeVhdAvailable) { throw 'Optimize-VHD is unavailable.' }
$principal = [Security.Principal.WindowsPrincipal][Security.Principal.WindowsIdentity]::GetCurrent()
if (-not $principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)) {
    throw 'Run the approved offline operation from Administrator Windows PowerShell.'
}
if (($file.Attributes -band [IO.FileAttributes]::SparseFile) -ne 0) {
    throw 'The registered VHDX is sparse; no storage settings will be changed automatically.'
}
if ($running) {
    if (-not $StopDistribution) { throw 'Save other Ubuntu work before explicitly approving StopDistribution.' }
    & wsl.exe --distribution $DistributionName --exec /bin/sync
    if ($LASTEXITCODE -ne 0) { throw 'Filesystem sync failed; distribution left running.' }
    & wsl.exe --terminate $DistributionName
    if ($LASTEXITCODE -ne 0) { throw 'Distribution did not stop.' }
    $deadline = [DateTime]::UtcNow.AddSeconds(30)
    while ($DistributionName -in (Get-RunningDistributions)) {
        if ([DateTime]::UtcNow -ge $deadline) { throw 'Distribution remains active; compaction refused.' }
        Start-Sleep -Milliseconds 250
    }
}
# The supported Hyper-V cmdlet requires the selected VHDX to be offline.
# It compacts unused blocks; it does not remove files or rewrite Git history.
Optimize-VHD -Path $file.FullName -Mode Full -ErrorAction Stop
$after = Get-Item -LiteralPath $file.FullName -Force
$report.LengthAfter = $after.Length
$report.HostFreeAfter = (Get-PSDrive -Name $drive).Free
$report.Completed = $true
$report | ConvertTo-Json
