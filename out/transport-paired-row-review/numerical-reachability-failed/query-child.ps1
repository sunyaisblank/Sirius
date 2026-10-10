$ErrorActionPreference='Stop'
[Console]::OutputEncoding=[System.Text.UTF8Encoding]::new($false)
$env:TEMP='\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\numerical\native-temp';$env:TMP=$env:TEMP
try {
 Add-Type -TypeDefinition ([IO.File]::ReadAllText('\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\numerical\native_gate.cs'))
 $sequence=Get-Content -Raw -LiteralPath '\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\numerical\sequences.json' | ConvertFrom-Json
 $result=[NativePairedTransportGate]::Run('C:\Windows\System32\vulkan-1.dll','\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\numerical\baseline.spv','\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\numerical\initial.bin','', '\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\numerical\observations\output',24,112776,'b52e6373f6b30aac6c7fe877f6e07f3c47c7ac9eafa1671840eacd19960d7cb1','316a6c591ba4e0b7087813ad91596217d5e3b76c2f5b4cf0e032a7ada4cf3b2c',$false,'\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\numerical\candidate.spv','721be7625f8eebbf722adf3e79a89050bffe9bf9b7260b8f5120b0a85ccd7f27',[string[]]@($sequence | ForEach-Object {$_.input}),[string[]]@($sequence | ForEach-Object {$_.input_sha256}),[System.UInt32[]]@($sequence | ForEach-Object {$_.active_rows}),[int[]]@($sequence | ForEach-Object {$_.upload_bytes}),[string[]]@($sequence | ForEach-Object {$_.name}))
 [IO.File]::WriteAllText('\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\numerical\query.json',($result | ConvertTo-Json -Depth 14),[Text.UTF8Encoding]::new($false))
 ConvertTo-Json -InputObject @($result['loaded_modules']) -Depth 8 | Set-Content -LiteralPath '\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\numerical\run-query\driver-modules.json' -Encoding UTF8
 [Environment]::Exit(0)
} catch { [Console]::Error.WriteLine($_.Exception.ToString());[Environment]::Exit(1) }
