$ErrorActionPreference='Stop'
[Console]::OutputEncoding=[System.Text.UTF8Encoding]::new($false)
$env:TEMP='\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\native-ordinary-shader-info\native-temp';$env:TMP=$env:TEMP
try {
 Add-Type -TypeDefinition ([IO.File]::ReadAllText('\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\native-ordinary-shader-info\native_census.cs'))
 $result=[NativeTransportFmaProbe]::Census('C:\Windows\System32\vulkan-1.dll')
 [IO.File]::WriteAllText('\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\native-ordinary-shader-info\census.json',($result | ConvertTo-Json -Depth 10),[Text.UTF8Encoding]::new($false))
 ConvertTo-Json -InputObject @($result['loaded_modules']) -Depth 8 | Set-Content -LiteralPath '\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\native-ordinary-shader-info\run\driver-modules.json' -Encoding UTF8
 [Console]::WriteLine('AMD shader-info advertised: '+$result['VK_AMD_shader_info_advertised'])
 [Environment]::Exit(0)
} catch { [Console]::Error.WriteLine($_.Exception.ToString());[Environment]::Exit(1) }
