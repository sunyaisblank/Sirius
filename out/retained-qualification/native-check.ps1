[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new($false)
$env:TEMP = '\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\retained-qualification\native-temp'
$env:TMP = $env:TEMP
$env:TMPDIR = '\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\renders\production-acceptance-b789ee3\native-radeon'
$env:SIRIUS_MEMORY_BUDGET_MB = '2048'
$env:SIRIUS_VULKAN_DEVICE = '0'
& '\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\bin\windows-msvc\production-acceptance-b789ee3\sirius_render_tests.exe' '--gtest_filter=VulkanRenderSession.RetainedMovingThinLensKerrDetectorMatchesCpuLinearRadiance' '--gtest_output=xml:\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\attestations\production-acceptance\b789ee3\radeon-moving-thinlens-timing.xml' '--gtest_repeat=1'
exit $LASTEXITCODE
