# Compile-only, revision-bound matched producer. No PE/provider/test execution.
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$PSNativeCommandUseErrorActionPreference = $false
$Root = [IO.Path]::GetFullPath($env:GITHUB_WORKSPACE)
Set-Location -LiteralPath $Root
$Build = Join-Path $Root 'bin/windows-msvc'
$Export = Join-Path $Build 'private-native-coefficient-sharing-export'
if (Test-Path -LiteralPath $Export) { throw 'Refuse an existing export; no stale products may be reused' }
New-Item -ItemType Directory -Path $Export | Out-Null
$Commands = [Collections.Generic.List[object]]::new()
$AllowedToolAliases = @{}
$Result = [ordered]@{
    schema = 'sirius-native-packed-matched-ci-producer-v1'
    producer_profile = 'native-coefficient-sharing-one-shader-only'
    canonical_CXX_source_revision = '2c8ef08693d46f57d5dea14efdbddd4874dac684'
    expected_private_compile_objects = 5
    passed = $false; configure_passed = $false; build_passed = $false
    source_unchanged = $false; frozen_build_inputs_unchanged = $false
    repository = $env:GITHUB_REPOSITORY; run_id = $env:GITHUB_RUN_ID
    run_attempt = $env:GITHUB_RUN_ATTEMPT; github_sha = $env:GITHUB_SHA
    source_revision = $null; commands = $Commands
    source_root = $Root; build_root = $Build; export_root = $Export
    pe_or_provider_or_test_executed = $false; mandatory_gate_executed = $false
    native_build_attestation_claimed = $false; scientific_qualification_claimed = $false
    release_or_adoption_claimed = $false; failure = $null
    scope = 'Compile-only two-arm matched products, with actual canonical dependency work. Input scope is tracked/private sources, fetched source trees, selected primary tools and adjacent DLLs; not complete host/SDK-header/DSO closure. Runtime/PE selectors and source/header/map proof require ROOT download review.'
}
function Write-Json($Value, [string]$Relative) {
    $Path = Join-Path $Export $Relative
    New-Item -ItemType Directory -Path (Split-Path -Parent $Path) -Force | Out-Null
    [IO.File]::WriteAllText($Path, ($Value | ConvertTo-Json -Depth 100) + "`n", [Text.UTF8Encoding]::new($false))
}
function File-Pin([string]$Path) {
    $Item = Get-Item -LiteralPath $Path -Force
    if ($Item.PSIsContainer) { throw "Expected file: $Path" }
    $Alias = ($Item.Attributes -band [IO.FileAttributes]::ReparsePoint) -or $Item.Length -eq 0
    if (($Item.Attributes -band [IO.FileAttributes]::ReparsePoint) -and -not $AllowedToolAliases.ContainsKey($Item.FullName)) {
        throw "Unadmitted reparse input/product: $Path"
    }
    $Stream = [IO.File]::OpenRead($Item.FullName)
    try { $ReadableBytes = $Stream.Length } finally { $Stream.Dispose() }
    return [ordered]@{ path = $Item.FullName; bytes = $Item.Length; readable_stream_bytes = $ReadableBytes
        sha256 = (Get-FileHash -LiteralPath $Item.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
        windows_attributes = $Item.Attributes.ToString(); link_type = [string]$Item.LinkType; link_target = @($Item.Target)
        configured_alias_metadata = [bool]($Alias -and $AllowedToolAliases.ContainsKey($Item.FullName))
        mode_scope = 'FileInfo length/attributes and readable stream hash recorded separately; alias metadata is not a backed interpreter or substitution authority' }
}
function Snapshot($Paths) {
    $Rows = [ordered]@{}
    foreach ($Path in @($Paths | Sort-Object -Unique)) { $Rows[$Path] = File-Pin $Path }
    return $Rows
}
function Assert-Same($Before, $After, [string]$Label) {
    if (($Before | ConvertTo-Json -Depth 100 -Compress) -cne ($After | ConvertTo-Json -Depth 100 -Compress)) {
        throw "$Label changed in membership, bytes or Windows attributes"
    }
}
function Invoke-Recorded([string]$Label, [string]$Tool, [string[]]$Arguments) {
    $Started = [DateTime]::UtcNow.ToString('o')
    & $Tool @Arguments 1> (Join-Path $Export "$Label.stdout.log") 2> (Join-Path $Export "$Label.stderr.log")
    $Code = $LASTEXITCODE
    $Record = [ordered]@{ label = $Label; argv = @($Tool) + @($Arguments); cwd = $Root
        started_utc = $Started; finished_utc = [DateTime]::UtcNow.ToString('o'); exit_code = $Code }
    $Commands.Add($Record)
    Write-Json $Commands 'commands.json'
    return $Code
}
function Clean-Revision {
    $Revision = (& $Git rev-parse HEAD).Trim()
    if ($LASTEXITCODE -ne 0 -or $Revision -cne $env:GITHUB_SHA) { throw 'HEAD differs from GITHUB_SHA' }
    $Status = @(& $Git status --porcelain --untracked-files=normal)
    if ($LASTEXITCODE -ne 0 -or $Status.Count -ne 0) { throw 'Git source is not clean' }
    return $Revision
}
function Source-Paths {
    $Names = ((@(& $Git -c core.quotepath=false ls-files -z) -join "`n") -split "`0") | Where-Object { $_ -ne '' }
    if ($LASTEXITCODE -ne 0) { throw 'Cannot enumerate tracked source' }
    $Paths = @($Names | ForEach-Object { Join-Path $Root $_ })
    foreach ($Relative in @('out/native-coefficient-sharing-producer-inputs',
                           'out/native-coefficient-sharing-windows-source')) {
        $Directory = Join-Path $Root $Relative
        if (-not (Test-Path -LiteralPath $Directory -PathType Container)) { throw "Missing private input tree: $Relative" }
        $Paths += @(Get-ChildItem -LiteralPath $Directory -Recurse -File -Force | ForEach-Object FullName)
    }
    return @($Paths | Sort-Object -Unique)
}
function Cache-Value([string]$Name) {
    $Found = [regex]::Matches($Cache, '(?m)^' + [regex]::Escape($Name) + ':[^=]+=([^\r\n]+)\r?$')
    $Values = @($Found | ForEach-Object { $_.Groups[1].Value } | Sort-Object -Unique)
    if ($Values.Count -ne 1) { throw "Absent or ambiguous configured tool/value: $Name" }
    return $Values[0]
}
function Copy-Evidence([string]$Path, [string]$Relative) {
    $Destination = Join-Path $Export $Relative
    New-Item -ItemType Directory -Path (Split-Path -Parent $Destination) -Force | Out-Null
    Copy-Item -LiteralPath $Path -Destination $Destination
    if ((File-Pin $Path).sha256 -cne (File-Pin $Destination).sha256) { throw "Copy mismatch: $Path" }
}
try {
    $GitCandidates = @(Get-Command git.exe -CommandType Application -ErrorAction Stop)
    $CMakeCandidates = @(Get-Command cmake.exe -CommandType Application -ErrorAction Stop)
    if ($GitCandidates.Count -eq 0 -or $CMakeCandidates.Count -eq 0) { throw 'Primary command resolution is empty' }
    $Git = [string]($GitCandidates[0].Source)
    $CMake = [string]($CMakeCandidates[0].Source)
    Write-Json ([ordered]@{
        git_candidates = @($GitCandidates | ForEach-Object Source); git_selected = $Git
        cmake_candidates = @($CMakeCandidates | ForEach-Object Source); cmake_selected = $CMake
        selection = 'First normal PATH Application entry as an explicit scalar string; PATH unchanged'
    }) 'command-resolution.json'
    $CTest = Join-Path (Split-Path -Parent $CMake) 'ctest.exe'
    $Fragment = Join-Path $Root 'out/native-coefficient-sharing-windows-source/matched-producer.cmake'
    $Result.source_revision = Clean-Revision
    if ($env:GITHUB_REF -cne 'refs/heads/development/native-coefficient-sharing-diagnostic-2026-10-20') { throw 'Wrong private diagnostic branch' }
    $SourceBeforeConfigure = Snapshot (Source-Paths)
    Write-Json $SourceBeforeConfigure 'source-before-configure.json'
    $Configure = @('--preset', 'windows-msvc', '-DSIRIUS_ALIGNMENT_MODE=qualification',
        '-DSIRIUS_REQUIRE_VULKAN_RUNTIME=ON', '-DCMAKE_CONFIGURATION_TYPES:STRING=Release',
        "-DCMAKE_PROJECT_Sirius_INCLUDE=$Fragment")
    if ((Invoke-Recorded 'configure' $CMake $Configure) -ne 0) { throw 'Configure failed; build not started' }
    $Result.configure_passed = $true
    $Cache = Get-Content -LiteralPath (Join-Path $Build 'CMakeCache.txt') -Raw
    # Preserve actual configure products before strict post-configure checks.
    Copy-Evidence (Join-Path $Build 'CMakeCache.txt') 'metadata/config/CMakeCache.txt'
    Copy-Evidence (Join-Path $Build 'CMakeFiles/CMakeConfigureLog.yaml') 'metadata/config/CMakeFiles/CMakeConfigureLog.yaml'
    foreach ($Required in @('BUILD_TESTS:BOOL=ON', 'SIRIUS_MANDATORY_TESTS:BOOL=ON',
                            'SIRIUS_WERROR:BOOL=ON', 'SIRIUS_CONTRACT_MODE:STRING=2',
                            'SIRIUS_REQUIRE_VULKAN_RUNTIME:BOOL=ON', 'SIRIUS_ALIGNMENT_MODE:STRING=qualification',
                            'CMAKE_CONFIGURATION_TYPES:STRING=Release', 'BUILD_SHARED_LIBS:BOOL=OFF')) {
        if ($Cache -notmatch ('(?m)^' + [regex]::Escape($Required) + '\r?$')) { throw "Strict cache mismatch: $Required" }
    }
    if ($Cache -match '(?m)^GLFW_LIBRARY_TYPE:STRING=(?!STATIC\r?$)[^\r\n]+') { throw 'Unrepresented GLFW DLL topology' }
    $CompilerFiles = @(Get-ChildItem -Path (Join-Path $Build 'CMakeFiles/*/CMakeCXXCompiler.cmake') -File)
    if ($CompilerFiles.Count -ne 1) { throw 'Ambiguous configured C++ compiler' }
    $CompilerFile = $CompilerFiles[0]
    $CompilerText = Get-Content -LiteralPath $CompilerFile.FullName -Raw
    $CompilerMatch = [regex]::Match($CompilerText, 'set\(CMAKE_CXX_COMPILER "([^"\r\n]+)"\)')
    $CompilerVersion = [regex]::Match($CompilerText, 'set\(CMAKE_CXX_COMPILER_VERSION "(19\.(\d+)\.[^"]+)"\)')
    if (-not $CompilerMatch.Success -or -not $CompilerVersion.Success) { throw 'Selected MSVC path/version absent' }
    $PythonPaths = @([regex]::Matches($Cache, '(?m)^(?:_Python3_EXECUTABLE|Python3_EXECUTABLE):[^=]+=([^\r\n]+)\r?$') | ForEach-Object { $_.Groups[1].Value } | Sort-Object -Unique)
    if ($PythonPaths.Count -ne 1) { throw 'Selected Python interpreter absent or ambiguous' }
    $ConfigureLog = Get-Content -LiteralPath (Join-Path $Build 'CMakeFiles/CMakeConfigureLog.yaml') -Raw
    $MSBuildPaths = @([regex]::Matches($ConfigureLog, 'Run Build Command\(s\):\s+"([^"\r\n]*[/\\]MSBuild\.exe)"') | ForEach-Object { $_.Groups[1].Value } | Sort-Object -Unique)
    if ($MSBuildPaths.Count -ne 1) { throw 'Actual configure MSBuild command absent or ambiguous; no guessed tool path allowed' }
    $PythonItem = Get-Item -LiteralPath $PythonPaths[0]
    $AllowedToolAliases[$PythonItem.FullName] = $true
    $Tools = [ordered]@{ cmake = $CMake; ctest = $CTest; python = $PythonPaths[0]; git = $Git
        powershell = (Join-Path $PSHOME 'pwsh.exe'); cl = $CompilerMatch.Groups[1].Value
        link = (Cache-Value 'CMAKE_LINKER'); lib = (Cache-Value 'CMAKE_AR')
        msbuild = $MSBuildPaths[0] }
    foreach ($Name in @('SIRIUS_SLANGC','SIRIUS_SPIRV_AS','SIRIUS_SPIRV_DIS','SIRIUS_SPIRV_VAL','SIRIUS_SPIRV_OPT')) {
        $Tools[$Name] = Cache-Value $Name
    }
    $PythonBinary = Join-Path (Split-Path -Parent $PythonPaths[0]) 'python.exe'
    if (-not (Test-Path -LiteralPath $PythonBinary -PathType Leaf) -or (Get-Item -LiteralPath $PythonBinary).Length -le 0) { throw 'Actual nonempty Python.exe beside configured interpreter is absent' }
    $Tools['python_nonempty_sibling'] = $PythonBinary
    Write-Json ([ordered]@{ configured_path = $PythonPaths[0]; nonempty_sibling = $PythonBinary
        alias_scope = 'Preserve actual configured alias length/reparse/target/stream data plus sibling bytes; no execution substitution or complete interpreter backup authority' }) 'python-alias-scope.json'
    $ToolPaths = @()
    foreach ($Tool in $Tools.Values) {
        $Selected = Get-Item -LiteralPath $Tool -ErrorAction Stop
        if ($Selected.Length -eq 0 -and -not $AllowedToolAliases.ContainsKey($Selected.FullName)) { throw "Empty primary tool: $Tool" }
        $ToolPaths += $Selected.FullName
        $ToolPaths += @(Get-ChildItem -LiteralPath (Split-Path -Parent $Tool) -Filter '*.dll' -File -Force | ForEach-Object FullName)
    }
    Write-Json $Tools 'selected-tools.json'
    $ConfiguredBuildInputs = @((Join-Path $Build 'CMakeCache.txt'), $CompilerFile.FullName,
        (Join-Path $CompilerFile.Directory.FullName 'CMakeCCompiler.cmake'),
        (Join-Path $CompilerFile.Directory.FullName 'CMakeSystem.cmake'),
        (Join-Path $Build 'CMakeFiles/CMakeConfigureLog.yaml'))
    $ConfiguredBuildInputs += @(Get-ChildItem -LiteralPath $Build -Recurse -File -Force | Where-Object { $_.Extension -in @('.vcxproj','.props','.targets','.sln') } | ForEach-Object FullName)
    function Build-Input-Paths {
        $Paths = @(Source-Paths) + @($ToolPaths) + @($ConfiguredBuildInputs)
        $Dependencies = Join-Path $Build '_deps'
        if (-not (Test-Path -LiteralPath $Dependencies -PathType Container)) { throw 'Fetched source directory absent' }
        $Trees = @(Get-ChildItem -LiteralPath $Dependencies -Directory -Filter '*-src' -Force)
        if ($Trees.Count -eq 0) { throw 'No fetched source trees' }
        foreach ($Tree in $Trees) { $Paths += @(Get-ChildItem -LiteralPath $Tree.FullName -Recurse -File -Force | ForEach-Object FullName) }
        return @($Paths | Sort-Object -Unique)
    }
    Assert-Same $SourceBeforeConfigure (Snapshot (Source-Paths)) 'Configure source/private inputs'
    $Before = Snapshot (Build-Input-Paths)
    Write-Json $Before 'build-inputs-before.json'
    foreach ($Gate in @('mandatory_gate.json','native_build_gate.json')) {
        if (Test-Path -LiteralPath (Join-Path $Build "generated/sirius/$Gate")) { throw "Unexpected promotion receipt: $Gate" }
    }
    $BuildArgs = @('--build', '--preset', 'windows-msvc', '--config', 'Release',
        '--target', 'SiriusNativeCoefficientSharingMatchedProducer', '--parallel', '2', '--',
        '/verbosity:normal', "/bl:$(Join-Path $Export 'build.binlog')")
    $BuildCode = Invoke-Recorded 'build' $CMake $BuildArgs
    $After = Snapshot (Build-Input-Paths)
    Write-Json $After 'build-inputs-after.json'
    Assert-Same $Before $After 'Frozen build inputs'
    $Result.frozen_build_inputs_unchanged = $true
    if ((Clean-Revision) -cne $Result.source_revision) { throw 'Build changed source revision' }
    Assert-Same $SourceBeforeConfigure (Snapshot (Source-Paths)) 'Complete source/private inputs'
    $Result.source_unchanged = $true
    $Result.build_passed = ($BuildCode -eq 0)
    foreach ($Gate in @('mandatory_gate.json','native_build_gate.json')) {
        if (Test-Path -LiteralPath (Join-Path $Build "generated/sirius/$Gate")) { throw "Private target produced forbidden promotion receipt: $Gate" }
    }
    # Preserve actual configured evidence/tlogs even when the build failed.
    $Config = @('CMakeCache.txt', 'CMakeFiles/CMakeConfigureLog.yaml', 'generated/sirius/alignment_receipt.json',
        'generated/sirius/base/operating_model_embedded.h')
    $VersionDirectory = [IO.Path]::GetRelativePath($Build, $CompilerFile.Directory.FullName)
    $Config += @("$VersionDirectory/CMakeCXXCompiler.cmake", "$VersionDirectory/CMakeCCompiler.cmake", "$VersionDirectory/CMakeSystem.cmake")
    foreach ($Relative in $Config) { Copy-Evidence (Join-Path $Build $Relative) "metadata/config/$Relative" }
    $Targets = @('native_coefficient_sharing_dopri','native_coefficient_sharing_baseline_compute','native_coefficient_sharing_candidate_compute',
        'native_coefficient_sharing_baseline_controls','native_coefficient_sharing_candidate_controls')
    foreach ($Target in $Targets) {
        Copy-Evidence (Join-Path $Build "$Target.vcxproj") "metadata/projects/$Target.vcxproj"
        $ObjectRoot = Join-Path $Build "$Target.dir/Release"
        if (Test-Path -LiteralPath $ObjectRoot -PathType Container) {
            foreach ($File in @(Get-ChildItem -LiteralPath $ObjectRoot -Recurse -File -Force | Where-Object { $_.Extension -in @('.tlog','.obj') })) {
                $Relative = [IO.Path]::GetRelativePath($Build, $File.FullName)
                Copy-Evidence $File.FullName "metadata/private-objects/$Relative"
            }
        }
    }
    # Export existing canonical compile/archive evidence; no extra build target.
    $BackendProject = 'src/sirius/backend/sirius_backend.vcxproj'
    Copy-Evidence (Join-Path $Build $BackendProject) "metadata/canonical-backend/$BackendProject"
    $BackendObjects = Join-Path $Build 'src/sirius/backend/sirius_backend.dir/Release'
    if (Test-Path -LiteralPath $BackendObjects -PathType Container) {
        foreach ($File in @(Get-ChildItem -LiteralPath $BackendObjects -Recurse -File -Force | Where-Object { $_.Extension -in @('.tlog','.obj') })) {
            $Relative = [IO.Path]::GetRelativePath($Build, $File.FullName)
            Copy-Evidence $File.FullName "metadata/canonical-backend/$Relative"
        }
    }
    $BackendArchive = 'src/sirius/backend/Release/sirius_backend.lib'
    if (Test-Path -LiteralPath (Join-Path $Build $BackendArchive) -PathType Leaf) {
        Copy-Evidence (Join-Path $Build $BackendArchive) "metadata/canonical-backend/$BackendArchive"
    }
    if ($BuildCode -ne 0) { throw 'Build failed; logs/input snapshots retained, no product PASS' }
    $Alignment = Get-Content -LiteralPath (Join-Path $Build 'generated/sirius/alignment_receipt.json') -Raw | ConvertFrom-Json
    if ($Alignment.source_revision -cne $Result.source_revision) { throw 'Alignment revision mismatch' }
    $Instance = Cache-Value 'CMAKE_GENERATOR_INSTANCE'
    $RedistVersionFile = Join-Path $Instance 'VC/Auxiliary/Build/Microsoft.VCRedistVersion.default.txt'
    $RedistVersion = (Get-Content -LiteralPath $RedistVersionFile -Raw).Trim()
    if ($RedistVersion -notmatch ('^14\.' + $CompilerVersion.Groups[2].Value + '\.')) { throw 'Selected MSVC/release CRT minor versions differ' }
    $CRTDirectory = Join-Path $Instance "VC/Redist/MSVC/$RedistVersion/x64/Microsoft.VC143.CRT"
    $CRTNames = @('concrt140.dll','msvcp140.dll','msvcp140_1.dll','msvcp140_2.dll',
        'msvcp140_atomic_wait.dll','msvcp140_codecvt_ids.dll','vcruntime140.dll',
        'vcruntime140_1.dll','Microsoft.VC143.CRT.manifest')
    $CRT = @($CRTNames | ForEach-Object { $Path = Join-Path $CRTDirectory $_; if (Test-Path -LiteralPath $Path -PathType Leaf) { Get-Item -LiteralPath $Path } })
    foreach ($Name in @('msvcp140.dll','vcruntime140.dll','vcruntime140_1.dll')) {
        if ($Name -notin $CRT.Name) { throw "Selected release CRT missing: $Name" }
    }
    $Loader = Get-Item -LiteralPath (Join-Path $env:VULKAN_SDK 'runtime/x64/vulkan-1.dll')
    $Origins = @{}
    foreach ($DLL in @($Loader) + @($CRT)) { $Origins[$DLL.Name] = File-Pin $DLL.FullName }
    Write-Json ([ordered]@{ msvc_version = $CompilerVersion.Groups[1].Value; release_crt_version = $RedistVersion
        redist_version_file = (File-Pin $RedistVersionFile); runtime_origins = $Origins }) 'runtime-packaging.json'
    foreach ($Arm in @('baseline','candidate')) {
        $Output = Join-Path $Build "private-native-coefficient-sharing/$Arm/Release"
        foreach ($Extension in @('exe','map')) {
            $Path = Join-Path $Output "original_${Arm}_controls.$Extension"
            Copy-Evidence $Path "payload/$Arm/original_${Arm}_controls.$Extension"
        }
        foreach ($DLL in @($Loader) + @($CRT)) { Copy-Evidence $DLL.FullName "payload/$Arm/$($DLL.Name)" }
    }
    $Result.private_compile_object_count = @((Get-ChildItem -LiteralPath (Join-Path $Export 'metadata/private-objects') -Recurse -File -Filter '*.obj')).Count
    if ($Result.private_compile_object_count -ne 5) { throw 'Expected exactly five private compile objects' }
    $Result.input_rows = $Before.Count
    $Result.passed = $true
} catch {
    $Result.failure = $_.Exception.ToString()
    $Result.failure_details = [ordered]@{
        ScriptStackTrace = $_.ScriptStackTrace
        InvocationInfo = [ordered]@{ PositionMessage = $_.InvocationInfo.PositionMessage }
    }
    Write-Error -Message $Result.failure -ErrorAction Continue
} finally {
    $Result.finished_utc = [DateTime]::UtcNow.ToString('o')
    Write-Json $Result 'producer-result.json'
    $Files = [ordered]@{}
    foreach ($File in @(Get-ChildItem -LiteralPath $Export -Recurse -File -Force)) {
        if ($File.Name -eq 'manifest.json') { continue }
        $Relative = [IO.Path]::GetRelativePath($Export, $File.FullName).Replace('\','/')
        $Files[$Relative] = File-Pin $File.FullName
    }
    Write-Json ([ordered]@{ schema = 'sirius-native-packed-ci-export-manifest-v1'; source_revision = $Result.source_revision
        passed = $Result.passed; files = $Files; complete_host_or_sdk_or_tool_dso_backup = $false }) 'manifest.json'
}
if (-not $Result.passed) { exit 1 }
