param([Parameter(Mandatory=$true)][string]$Executable,[Parameter(Mandatory=$true)][string]$Output)
$ErrorActionPreference='Stop'
New-Item -ItemType Directory -Path $Output -Force | Out-Null
$identity=[ordered]@{platform='Windows';source_revision='6e7ba9347a6eb15ccc607b85917edb7d9be5c34b';qualification_claimed=$false;executable=$Executable;sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $Executable).Hash.ToLowerInvariant();captured_utc=[DateTime]::UtcNow.ToString('o');command=@('--json','info','system')}
& $Executable --json info system 1> (Join-Path $Output 'inventory.json') 2> (Join-Path $Output 'inventory-stderr.log')
$identity.returncode=$LASTEXITCODE
$identity | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'identity.json') -Encoding UTF8
if ($LASTEXITCODE -ne 0) { throw 'Actual native inventory failed' }
$loader=Join-Path $env:SystemRoot 'System32\vulkan-1.dll'
if (Test-Path $loader) { [ordered]@{path=$loader;sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $loader).Hash.ToLowerInvariant();version=(Get-Item -LiteralPath $loader).VersionInfo.FileVersion;scope='Installed loader input identity, not actual loaded-module proof'} | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $Output 'loader-input.json') -Encoding UTF8 }
Write-Output 'Actual native inventory captured'
