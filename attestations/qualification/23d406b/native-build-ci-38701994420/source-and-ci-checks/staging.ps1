$ErrorActionPreference = "Stop"
$revision = (git rev-parse HEAD).Trim()
if ($LASTEXITCODE -ne 0 -or $revision -ne $env:GITHUB_SHA -or (git status --porcelain)) { throw "diagnostic source identity mismatch" }
$alignment = "bin/windows-msvc/generated/sirius/alignment_receipt.json"
if ((Get-Content -Raw $alignment | ConvertFrom-Json).source_revision -ne $revision) { throw "diagnostic alignment revision mismatch" }
$volumes = @("bin/windows-msvc/src/sirius/app/Release", "bin/windows-msvc/tests/backend/Release")
$executables = @("sirius.exe", "sirius_backend_tests.exe")
$loader = Get-Item -LiteralPath (Join-Path $env:VULKAN_SDK "runtime/x64/vulkan-1.dll")
$cache = Get-Content -Raw "bin/windows-msvc/CMakeCache.txt"
if ($cache -notmatch '(?m)^BUILD_SHARED_LIBS:BOOL=OFF\r?$') { throw "diagnostic volumes require the existing static-library topology" }
if ($cache -match '(?m)^GLFW_LIBRARY_TYPE:STRING=(?!STATIC\r?$)[^\r\n]+') { throw "diagnostic GLFW DLL topology is not represented" }
$instance = [regex]::Match($cache, '(?m)^CMAKE_GENERATOR_INSTANCE:[^=]+=(.+)$')
if (-not $instance.Success) { throw "selected Visual Studio instance is absent" }
$compiler = @(Get-ChildItem -Path "bin/windows-msvc/CMakeFiles/*/CMakeCXXCompiler.cmake" -File)
if ($compiler.Count -ne 1) { throw "selected C++ compiler metadata is ambiguous" }
$compilerVersion = [regex]::Match((Get-Content -Raw $compiler[0].FullName), 'set\(CMAKE_CXX_COMPILER_VERSION "(19\.(\d+)\.[^"]+)"\)')
if (-not $compilerVersion.Success) { throw "selected MSVC compiler version is absent" }
$redistVersionFile = Join-Path ($instance.Groups[1].Value.Trim()) "VC/Auxiliary/Build/Microsoft.VCRedistVersion.default.txt"
$redistVersion = (Get-Content -Raw $redistVersionFile).Trim()
if ($redistVersion -notmatch ('^14\.' + $compilerVersion.Groups[2].Value + '\.')) { throw "selected MSVC and release CRT versions do not match" }
$crtDirectory = Join-Path ($instance.Groups[1].Value.Trim()) "VC/Redist/MSVC/$redistVersion/x64/Microsoft.VC143.CRT"
$crtNames = @("concrt140.dll", "msvcp140.dll", "msvcp140_1.dll", "msvcp140_2.dll", "msvcp140_atomic_wait.dll", "msvcp140_codecvt_ids.dll", "vcruntime140.dll", "vcruntime140_1.dll", "Microsoft.VC143.CRT.manifest")
$crt = @($crtNames | ForEach-Object { $path = Join-Path $crtDirectory $_; if (Test-Path -LiteralPath $path -PathType Leaf) { Get-Item -LiteralPath $path } })
foreach ($name in @("msvcp140.dll", "vcruntime140.dll", "vcruntime140_1.dll")) {
  if ($name -notin $crt.Name) { throw "selected release CRT missing: $name" }
}
$files = @{}
$origins = @{}
foreach ($dll in @($loader) + @($crt)) { $origins[$dll.Name] = $dll.FullName }
for ($index = 0; $index -lt $volumes.Count; $index++) {
  $volume = $volumes[$index]
  foreach ($dll in @($loader) + @($crt)) { Copy-Item -LiteralPath $dll.FullName -Destination $volume -Force }
  $payload = @(Get-Item -LiteralPath (Join-Path $volume $executables[$index]))
  $payload += @(Get-ChildItem -Path "$volume/*.dll", "$volume/*.manifest" -File -Force)
  $payload += @(Get-ChildItem -LiteralPath "$volume/resources" -Recurse -File -Force)
  foreach ($file in $payload) {
    $relative = [IO.Path]::GetRelativePath($PWD.Path, $file.FullName).Replace('\', '/')
    $files[$relative] = @{ bytes = $file.Length; sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $file.FullName).Hash.ToLowerInvariant() }
  }
}
$compilerCopy = "bin/windows-msvc/generated/sirius/windows-diagnostic-CMakeCXXCompiler.cmake"
Copy-Item -LiteralPath $compiler[0].FullName -Destination $compilerCopy
foreach ($path in @($alignment, "bin/windows-msvc/CMakeCache.txt", $compilerCopy)) {
  $file = Get-Item -LiteralPath $path
  $files[$path] = @{ bytes = $file.Length; sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $path).Hash.ToLowerInvariant() }
}
$record = @{ schema = "sirius-windows-diagnostic-build-products-v1"; source_revision = $revision; repository = $env:GITHUB_REPOSITORY; run_id = $env:GITHUB_RUN_ID; run_attempt = $env:GITHUB_RUN_ATTEMPT; msvc_version = $compilerVersion.Groups[1].Value; release_crt_version = $redistVersion; runtime_dll_origins = $origins; files = $files; diagnostic_build_only = $true; qualification_claimed = $false; scope = "Strict compiled CLI/backend-test volumes only; no native runtime, rendering, UI, scientific, release or native-build attestation claim. Windows OS/UCRT and the physical graphics driver remain host dependencies." }
$output = "bin/windows-msvc/generated/sirius/windows-diagnostic-build-products.json"
[IO.File]::WriteAllText((Join-Path $PWD.Path $output), ($record | ConvertTo-Json -Depth 8) + "`n", [Text.UTF8Encoding]::new($false))
