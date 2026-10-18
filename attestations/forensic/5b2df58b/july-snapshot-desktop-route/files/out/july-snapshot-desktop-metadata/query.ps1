# SOURCE ONLY. Fixed snapshot metadata paths; ROOT owns the external 60s bound.
$ErrorActionPreference = 'Stop'
[Console]::OutputEncoding = New-Object System.Text.UTF8Encoding($false)
$watch = [Diagnostics.Stopwatch]::StartNew()
$process = [Diagnostics.Process]::GetCurrentProcess()
$heldProcessHandle = $process.Handle
$script:entryCount = 0
$script:outputCount = 0
$script:callCount = 0
$script:api = $null
$script:findData = [IntPtr]::Zero
$script:cleanupErrors = New-Object 'Collections.Generic.List[string]'
$snapshots = @('HarddiskVolumeShadowCopy1', 'HarddiskVolumeShadowCopy2', 'HarddiskVolumeShadowCopy3')
$metadataPattern = '(?i)july|tape|manifest|comparator|forensic|sirius|dngr|readme|catalog|provenance|checksum'
$result = [ordered]@{
    schema = 'sirius-july-snapshot-desktop-metadata-v1'
    utc = [DateTime]::UtcNow.ToString('o')
    owner_pid = $PID
    owner_start_utc = $process.StartTime.ToUniversalTime().ToString('o')
    held_own_process_handle = $true
    soft_seconds = 30
    outer_seconds_required = 60
    maximum_entries = 4096
    maximum_output_entries = 256
    snapshots = @()
    gaps = @()
    read_only_metadata = $true
    scope = 'Only three observed snapshot targets: prefix attributes, Desktop one level, then direct Sirius/Remediation one level if present. No recursion, file contents, images, restore, mount, elevation, privilege changes, or original July identity conclusion.'
}
$base = Split-Path -Parent $MyInvocation.MyCommand.Path
@{schema='sirius-july-snapshot-desktop-bootstrap-v1';owner_pid=$PID;owner_start_utc=$result.owner_start_utc;held_own_process_handle=$true;soft_seconds=30;outer_seconds_required=60} | ConvertTo-Json -Compress | Set-Content -LiteralPath (Join-Path $base 'query-bootstrap.json') -Encoding UTF8
function Check-Bound {
    if ($watch.Elapsed.TotalSeconds -ge 30) { throw 'Soft 30-second limit reached between native calls' }
    if ($script:entryCount -ge 4096) { throw '4096 total-entry limit reached' }
}
function U32([int]$Offset) {
    return [BitConverter]::ToUInt32([BitConverter]::GetBytes([Runtime.InteropServices.Marshal]::ReadInt32($script:findData, $Offset)), 0)
}
function Current-Entry {
    # WIN32_FIND_DATAW is 592 bytes on the required x64 Windows ABI.
    $nameField = [Runtime.InteropServices.Marshal]::PtrToStringUni([IntPtr]($script:findData.ToInt64() + 44), 260)
    $terminator = $nameField.IndexOf([char]0)
    if ($terminator -lt 0) { throw 'Native filename is not terminated inside WIN32_FIND_DATAW' }
    $name = $nameField.Substring(0, $terminator)
    $attributes = U32 0
    $filetime = [Runtime.InteropServices.Marshal]::ReadInt64($script:findData, 20)
    $writeUtc = $null
    try { $writeUtc = [DateTime]::FromFileTimeUtc($filetime).ToString('o') } catch {}
    $script:entryCount++
    return @{
        name = $name
        attributes_hex = ('{0:X8}' -f $attributes)
        is_directory = (($attributes -band 16) -ne 0)
        is_reparse_point = (($attributes -band 1024) -ne 0)
        bytes_decimal = (([decimal](U32 28) * 4294967296) + [decimal](U32 32)).ToString([Globalization.CultureInfo]::InvariantCulture)
        last_write_filetime_hex = ('{0:X16}' -f $filetime)
        last_write_utc = $writeUtc
    }
}
function Close-Find([IntPtr]$Handle) {
    if ($Handle -eq [IntPtr](-1) -or $Handle -eq [IntPtr]::Zero) { return }
    $ok = [int]$script:api.GetMethod('FindClose').Invoke($null, ([object[]]@($Handle)))
    $errorCode = [Runtime.InteropServices.Marshal]::GetLastWin32Error()
    if ($ok -eq 0) { $script:cleanupErrors.Add('FindClose Win32=' + $errorCode) }
}
function Open-Find([string]$Path) {
    Check-Bound
    $pathBuffer = [Runtime.InteropServices.Marshal]::StringToHGlobalUni($Path)
    try {
        $handle = [IntPtr]$script:api.GetMethod('FindFirstFileW').Invoke($null, ([object[]]@($pathBuffer, $script:findData)))
        $errorCode = [Runtime.InteropServices.Marshal]::GetLastWin32Error()
        $script:callCount++
        return @{handle=$handle;error=$errorCode}
    } finally { [Runtime.InteropServices.Marshal]::FreeHGlobal($pathBuffer) }
}
function Exact-Directory([string]$Path) {
    $opened = Open-Find $Path
    if ($opened.handle -eq [IntPtr](-1)) {
        $status = if ($opened.error -eq 5) {'access_denied'} elseif ($opened.error -eq 2 -or $opened.error -eq 3) {'not_found'} else {'query_failed'}
        return @{path=$Path;status=$status;win32_error=$opened.error}
    }
    if ($opened.handle -eq [IntPtr]::Zero) { throw 'Unexpected null FindFirstFileW handle' }
    try {
        $entry = Current-Entry
        $status = if ($entry.is_reparse_point) {'reparse_not_followed'} elseif (-not $entry.is_directory) {'not_a_directory'} else {'directory_metadata_observed'}
        return @{path=$Path;status=$status;entry=$entry}
    } finally { Close-Find $opened.handle }
}
function List-One-Level([string]$Path, [bool]$Desktop) {
    $rows = New-Object 'Collections.Generic.List[object]'
    $directTargets = New-Object 'Collections.Generic.List[string]'
    $state = @{path=$Path;status='partial';enumeration_complete=$false;inspected=0;entries=@();direct_targets=@()}
    $opened = Open-Find ($Path + '\*')
    if ($opened.handle -eq [IntPtr](-1)) {
        $state.win32_error = $opened.error
        # FindFirst ERROR_FILE_NOT_FOUND on this already-validated directory can mean no wildcard entries.
        if ($opened.error -eq 2) { $state.status='empty';$state.enumeration_complete=$true }
        elseif ($opened.error -eq 5) { $state.status='access_denied' }
        else { $state.status='query_failed' }
        return $state
    }
    if ($opened.handle -eq [IntPtr]::Zero) { throw 'Unexpected null enumeration handle' }
    try {
        while ($true) {
            $entry = Current-Entry
            $state.inspected++
            if ($entry.name -ne '.' -and $entry.name -ne '..') {
                $knownDirectory = $entry.is_directory -and ($entry.name -ieq 'Sirius' -or $entry.name -ieq 'Remediation')
                $report = if ($Desktop) {$knownDirectory -or $entry.name -match $metadataPattern} else {$entry.is_directory -or $entry.name -match $metadataPattern}
                if ($report) {
                    if ($script:outputCount -ge 256) { throw '256 total-output-entry limit reached' }
                    $rows.Add($entry);$script:outputCount++
                }
                if ($Desktop -and $knownDirectory -and -not $entry.is_reparse_point) { $directTargets.Add($entry.name) }
            }
            Check-Bound
            $ok = [int]$script:api.GetMethod('FindNextFileW').Invoke($null, ([object[]]@($opened.handle, $script:findData)))
            $errorCode = [Runtime.InteropServices.Marshal]::GetLastWin32Error()
            $script:callCount++
            if ($ok -eq 0) {
                $state.win32_error = $errorCode
                if ($errorCode -eq 18) {$state.status='completed';$state.enumeration_complete=$true}
                elseif ($errorCode -eq 5) {$state.status='partial_access_denied'}
                else {$state.status='partial_query_failed'}
                break
            }
        }
    } catch { $state.gap=$_.Exception.Message }
    finally { Close-Find $opened.handle }
    $state.entries = $rows.ToArray()
    $state.direct_targets = $directTargets.ToArray()
    return $state
}
try {
    if ([IntPtr]::Size -ne 8 -or [Environment]::OSVersion.Platform -ne [PlatformID]::Win32NT -or $PSVersionTable.PSEdition -ne 'Desktop') { throw 'Requires existing x64 Windows PowerShell Desktop' }
    $assemblyName = New-Object Reflection.AssemblyName ('SiriusSnapshotDesktop' + $PID)
    $assembly = [AppDomain]::CurrentDomain.DefineDynamicAssembly($assemblyName, [Reflection.Emit.AssemblyBuilderAccess]::Run)
    $module = $assembly.DefineDynamicModule($assemblyName.Name)
    $builder = $module.DefineType('SiriusSnapshotDesktop', [Reflection.TypeAttributes]'Public, Abstract, Sealed')
    $dllImport = [Runtime.InteropServices.DllImportAttribute]
    $constructor = $dllImport.GetConstructor([Type[]]@([string]))
    $fields = [Reflection.FieldInfo[]]@($dllImport.GetField('EntryPoint'),$dllImport.GetField('ExactSpelling'),$dllImport.GetField('SetLastError'),$dllImport.GetField('CallingConvention'),$dllImport.GetField('CharSet'),$dllImport.GetField('PreserveSig'))
    foreach ($signature in @(
        @{Name='FindFirstFileW';Return=[IntPtr];Types=[Type[]]@([IntPtr],[IntPtr])},
        @{Name='FindNextFileW';Return=[int];Types=[Type[]]@([IntPtr],[IntPtr])},
        @{Name='FindClose';Return=[int];Types=[Type[]]@([IntPtr])}
    )) {
        $method = $builder.DefineMethod($signature.Name, [Reflection.MethodAttributes]'Public, Static, PinvokeImpl', [Reflection.CallingConventions]::Standard, $signature.Return, $signature.Types)
        $values = [object[]]@($signature.Name,$true,$true,[Runtime.InteropServices.CallingConvention]::Winapi,[Runtime.InteropServices.CharSet]::Unicode,$true)
        $attribute = [Reflection.Emit.CustomAttributeBuilder]::new($constructor,([object[]]@('kernel32.dll')),$fields,$values)
        $method.SetCustomAttribute($attribute)
        $method.SetImplementationFlags([Reflection.MethodImplAttributes]::PreserveSig)
    }
    $script:api = $builder.CreateType()
    $script:findData = [Runtime.InteropServices.Marshal]::AllocHGlobal(592)
    foreach ($snapshot in $snapshots) {
        Check-Bound
        $root = '\\?\GLOBALROOT\Device\' + $snapshot
        $item = @{snapshot=$snapshot;prefix_checks=@();desktop=$null;direct_directories=@();status='partial'}
        $result.snapshots += $item
        $prefixOkay = $true
        foreach ($relative in @('Users','Users\alarisadmin','Users\alarisadmin\Desktop')) {
            $checked = Exact-Directory ($root + '\' + $relative)
            $item.prefix_checks += $checked
            if ($checked.status -ne 'directory_metadata_observed') {$prefixOkay=$false;$item.status=$checked.status;break}
        }
        if (-not $prefixOkay) { continue }
        $desktopPath = $root + '\Users\alarisadmin\Desktop'
        $item.desktop = List-One-Level $desktopPath $true
        if (-not $item.desktop.enumeration_complete) { continue }
        foreach ($name in $item.desktop.direct_targets) {
            Check-Bound
            $path = $desktopPath + '\' + $name
            $checked = Exact-Directory $path
            if ($checked.status -ne 'directory_metadata_observed') {$item.direct_directories += @{path=$path;status=$checked.status;metadata=$checked};continue}
            $item.direct_directories += (List-One-Level $path $false)
        }
        $item.status = if (@($item.direct_directories | Where-Object {-not $_.enumeration_complete}).Count -eq 0) {'metadata_walk_finished'} else {'partial'}
    }
} catch { $result.gaps += $_.Exception.Message }
finally {
    if ($script:findData -ne [IntPtr]::Zero) {[Runtime.InteropServices.Marshal]::FreeHGlobal($script:findData)}
    $result.inspected_entries = $script:entryCount
    $result.output_entries = $script:outputCount
    $result.native_find_calls = $script:callCount
    $result.cleanup_errors = $script:cleanupErrors.ToArray()
    $result.elapsed_seconds = $watch.Elapsed.TotalSeconds
    $result.all_three_desktops_enumerated = ($result.snapshots.Count -eq 3 -and @($result.snapshots | Where-Object {-not $_.desktop -or -not $_.desktop.enumeration_complete}).Count -eq 0)
    $result.all_selected_existing_directories_enumerated = ($result.all_three_desktops_enumerated -and @($result.snapshots | Where-Object {$_.status -ne 'metadata_walk_finished'}).Count -eq 0)
    $result.audit_success = ($result.gaps.Count -eq 0 -and $script:cleanupErrors.Count -eq 0 -and $result.all_selected_existing_directories_enumerated)
    $result | ConvertTo-Json -Depth 10 -Compress
    $process.Dispose()
}
if (-not $result.audit_success) {exit 2}
