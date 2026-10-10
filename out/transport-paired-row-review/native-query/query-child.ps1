$ErrorActionPreference='Stop'
[Console]::OutputEncoding=[System.Text.UTF8Encoding]::new($false)
$env:TEMP='\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\native-query\query-native-temp';$env:TMP=$env:TEMP
try {
 Add-Type -TypeDefinition ([IO.File]::ReadAllText('\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\native-query\native_query.cs'))
 $result=[NativeTransportFmaProbe]::Observe('C:\Windows\System32\vulkan-1.dll','\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\paired-transport.spv','\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\paired-transport.spv','\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\native-query\payloads')
 [IO.File]::WriteAllText('\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\native-query\query.json',($result | ConvertTo-Json -Depth 14),[Text.UTF8Encoding]::new($false))
 ConvertTo-Json -InputObject @($result['loaded_modules']) -Depth 8 | Set-Content -LiteralPath '\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\transport-paired-row-review\native-query\run-query\driver-modules.json' -Encoding UTF8
 [Environment]::Exit(0)
} catch { [Console]::Error.WriteLine($_.Exception.ToString());[Environment]::Exit(1) }
