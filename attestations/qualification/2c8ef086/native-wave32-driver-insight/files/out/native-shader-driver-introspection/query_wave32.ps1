param([Parameter(Mandatory=$true)][string]$Config)
$ErrorActionPreference='Stop'
$settings=Get-Content -Raw -LiteralPath $Config | ConvertFrom-Json
foreach($name in @('VK_DRIVER_FILES','VK_ICD_FILENAMES','VK_ADD_DRIVER_FILES','VK_INSTANCE_LAYERS','AMD_VULKAN_ICD')) {
 if (![string]::IsNullOrEmpty([Environment]::GetEnvironmentVariable($name))) { throw ('Unexpected query override: '+$name) }
}
foreach($queryInput in @($settings.inputs)) {
 $item=Get-Item -LiteralPath $queryInput.path
 if ($item.Length -ne $queryInput.bytes -or (Get-FileHash -Algorithm SHA256 -LiteralPath $queryInput.path).Hash.ToLowerInvariant() -ne $queryInput.sha256) { throw 'Query input identity mismatch' }
}
Add-Type -Path $settings.csharp
$result=[SiriusShaderWave32Query]::Run($settings.loader,$settings.shader,$settings.query_output)
$result | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $settings.result -Encoding UTF8
foreach($queryInput in @($settings.inputs)) {
 if ((Get-FileHash -Algorithm SHA256 -LiteralPath $queryInput.path).Hash.ToLowerInvariant() -ne $queryInput.sha256) { throw 'Query input changed' }
}
