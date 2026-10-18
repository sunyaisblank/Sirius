# SOURCE ONLY: fixed, read-only NT object-directory name enumeration.
# ROOT must independently review and own a 60-second external process bound.
$ErrorActionPreference = 'Stop'
[Console]::OutputEncoding = New-Object System.Text.UTF8Encoding($false)
$watch = [Diagnostics.Stopwatch]::StartNew()
$process = [Diagnostics.Process]::GetCurrentProcess()
$heldProcessHandle = $process.Handle
$result = [ordered]@{
    schema = 'sirius-july-nt-snapshot-name-query-v1'
    utc = [DateTime]::UtcNow.ToString('o')
    owner_pid = $PID
    owner_start_utc = $process.StartTime.ToUniversalTime().ToString('o')
    held_own_process_handle = $true
    directory = '\Device'
    name_pattern = '^HarddiskVolumeShadowCopy[0-9]+$'
    requested_access = 'DIRECTORY_QUERY (0x00000001)'
    soft_seconds = 30
    outer_seconds_required = 60
    maximum_query_calls = 4096
    maximum_matching_names = 128
    native_buffer_bytes = 65536
    query_calls = 0
    inspected_entries = 0
    matches = @()
    open_status = $null
    last_query_status = $null
    close_status = $null
    end_of_directory_observed = $false
    directory_handle_closed = $false
    query_completed = $false
    coverage = 'unavailable'
    gaps = @()
    read_only_names = $true
    scope = 'One nonrecursive NT object-directory walk. No device/volume/file opening, snapshot contents, CIM/VSS query, restore, elevation, drive scan, or original July artifact identity claim. Enumeration is not atomic; zero observed matching names is not proof that no snapshots exist.'
}
$bootstrap = Join-Path (Split-Path -Parent $MyInvocation.MyCommand.Path) 'query-bootstrap.json'
@{schema='sirius-july-nt-snapshot-name-bootstrap-v1';owner_pid=$PID;owner_start_utc=$result.owner_start_utc;held_own_process_handle=$true;soft_seconds=30;outer_seconds_required=60} | ConvertTo-Json -Compress | Set-Content -LiteralPath $bootstrap -Encoding UTF8
$directoryHandle = [IntPtr]::Zero
$nameBuffer = [IntPtr]::Zero
$unicodeString = [IntPtr]::Zero
$attributes = [IntPtr]::Zero
$queryBuffer = [IntPtr]::Zero
$api = $null
$hits = New-Object 'Collections.Generic.List[object]'
function Status-Hex([int]$Status) { return ('{0:X8}' -f $Status) }
function Read-DirectoryString([IntPtr]$Buffer, [uint32]$ReturnedBytes, [int]$Offset) {
    # x64 UNICODE_STRING: USHORT Length, USHORT MaximumLength, padding, pointer.
    $length = [Runtime.InteropServices.Marshal]::ReadInt16($Buffer, $Offset) -band 65535
    $maximum = [Runtime.InteropServices.Marshal]::ReadInt16($Buffer, $Offset + 2) -band 65535
    $pointer = [Runtime.InteropServices.Marshal]::ReadIntPtr($Buffer, $Offset + 8)
    if (($length % 2) -ne 0 -or $maximum -lt $length) { throw 'Invalid returned UNICODE_STRING lengths' }
    if ($length -eq 0) { return '' }
    $first = $Buffer.ToInt64()
    $address = $pointer.ToInt64()
    if ($address -lt ($first + 32) -or ($address + $length) -gt ($first + $ReturnedBytes)) { throw 'Returned string lies outside the bounded native result buffer' }
    return [Runtime.InteropServices.Marshal]::PtrToStringUni($pointer, ($length / 2))
}
try {
    if ([IntPtr]::Size -ne 8 -or [Environment]::OSVersion.Platform -ne [PlatformID]::Win32NT) { throw 'Requires 64-bit Windows PowerShell' }
    if ($PSVersionTable.PSEdition -ne 'Desktop') { throw 'Requires existing Windows PowerShell .NET Framework Reflection.Emit API' }
    # Dynamic P/Invoke stubs only: no Add-Type, C# source, external csc, or child process.
    $assemblyName = New-Object Reflection.AssemblyName ('SiriusNtSnapshotNames' + $PID)
    $assembly = [AppDomain]::CurrentDomain.DefineDynamicAssembly($assemblyName, [Reflection.Emit.AssemblyBuilderAccess]::Run)
    $module = $assembly.DefineDynamicModule($assemblyName.Name)
    $builder = $module.DefineType('SiriusNtSnapshotNames', [Reflection.TypeAttributes]'Public, Abstract, Sealed')
    $signatures = @(
        @{Name='NtOpenDirectoryObject';Types=[Type[]]@([IntPtr].MakeByRefType(), [uint32], [IntPtr])},
        @{Name='NtQueryDirectoryObject';Types=[Type[]]@([IntPtr], [IntPtr], [uint32], [byte], [byte], [uint32].MakeByRefType(), [uint32].MakeByRefType())},
        @{Name='NtClose';Types=[Type[]]@([IntPtr])}
    )
    # All strings are manually built UTF-16 structures; P/Invoke receives only numeric/pointer arguments.
    # Ansi name resolution probes the exact exported Nt* spelling before any suffix.
    foreach ($signature in $signatures) {
        $method = $builder.DefinePInvokeMethod($signature.Name, 'ntdll.dll', $signature.Name, [Reflection.MethodAttributes]'Public, Static, PinvokeImpl', [Reflection.CallingConventions]::Standard, [int], $signature.Types, [Runtime.InteropServices.CallingConvention]::Winapi, [Runtime.InteropServices.CharSet]::Ansi)
        $method.SetImplementationFlags([Reflection.MethodImplAttributes]($method.GetMethodImplementationFlags() -bor [Reflection.MethodImplAttributes]::PreserveSig))
    }
    $api = $builder.CreateType()
    $nameBuffer = [Runtime.InteropServices.Marshal]::StringToHGlobalUni('\Device')
    $unicodeString = [Runtime.InteropServices.Marshal]::AllocHGlobal(16)
    [Runtime.InteropServices.Marshal]::Copy((New-Object byte[] 16), 0, $unicodeString, 16)
    [Runtime.InteropServices.Marshal]::WriteInt16($unicodeString, 0, 14)
    [Runtime.InteropServices.Marshal]::WriteInt16($unicodeString, 2, 16)
    [Runtime.InteropServices.Marshal]::WriteIntPtr($unicodeString, 8, $nameBuffer)
    $attributes = [Runtime.InteropServices.Marshal]::AllocHGlobal(48)
    [Runtime.InteropServices.Marshal]::Copy((New-Object byte[] 48), 0, $attributes, 48)
    [Runtime.InteropServices.Marshal]::WriteInt32($attributes, 0, 48)
    [Runtime.InteropServices.Marshal]::WriteIntPtr($attributes, 16, $unicodeString)
    [Runtime.InteropServices.Marshal]::WriteInt32($attributes, 24, 64) # OBJ_CASE_INSENSITIVE
    $openArgs = [object[]]@([IntPtr]::Zero, [uint32]1, $attributes)
    $status = [int]$api.GetMethod('NtOpenDirectoryObject').Invoke($null, $openArgs)
    $directoryHandle = [IntPtr]$openArgs[0]
    $result.open_status = Status-Hex $status
    if ($status -ne 0) {
        $result.coverage = if ($result.open_status -eq 'C0000022') { 'access_denied' } else { 'open_failed' }
        $result.gaps += 'NT directory query was not admitted; no snapshot absence conclusion.'
    } else {
        if ($directoryHandle -eq [IntPtr]::Zero) { throw 'Successful open returned a null handle' }
        $result.coverage = 'partial'
        $queryBuffer = [Runtime.InteropServices.Marshal]::AllocHGlobal(65536)
        $context = [uint32]0
        $restart = [byte]1
        while ($result.query_calls -lt 4096) {
            if ($watch.Elapsed.TotalSeconds -ge 30) { $result.gaps += 'Soft 30-second limit reached between native calls'; break }
            $queryArgs = [object[]]@($directoryHandle, $queryBuffer, [uint32]65536, [byte]1, $restart, $context, [uint32]0)
            $status = [int]$api.GetMethod('NtQueryDirectoryObject').Invoke($null, $queryArgs)
            $result.query_calls++
            $restart = [byte]0
            $context = [uint32]$queryArgs[5]
            $returned = [uint32]$queryArgs[6]
            $result.last_query_status = Status-Hex $status
            if ($result.last_query_status -eq '8000001A') { $result.end_of_directory_observed = $true; $result.query_completed = $true; $result.coverage = 'completed_visible_directory'; break }
            if ($status -ne 0 -and $result.last_query_status -ne '00000105') {
                $result.gaps += ('Native enumeration stopped: NTSTATUS ' + $result.last_query_status)
                if ($result.last_query_status -eq 'C0000022') { $result.coverage = 'partial_access_denied' }
                break
            }
            if ($returned -lt 32 -or $returned -gt 65536) { throw 'Native returned length outside the bounded OBJECT_DIRECTORY_INFORMATION result' }
            $entryName = Read-DirectoryString $queryBuffer $returned 0
            $entryType = Read-DirectoryString $queryBuffer $returned 16
            $result.inspected_entries++
            if ($entryName -imatch '^HarddiskVolumeShadowCopy[0-9]+$') {
                if ($hits.Count -ge 128) { $result.gaps += 'Matching-name output cap reached'; break }
                $hits.Add(@{name=$entryName;object_type=$entryType;nt_object_path=('\Device\' + $entryName)})
            }
        }
        if (-not $result.query_completed -and $result.query_calls -ge 4096) { $result.gaps += '4096 native-query call cap reached' }
    }
} catch {
    if ($result.coverage -eq 'unavailable') { $result.coverage = 'api_or_runtime_unavailable' }
    $result.gaps += $_.Exception.Message
} finally {
    if ($directoryHandle -ne [IntPtr]::Zero -and $api) {
        try {
            $closeStatus = [int]$api.GetMethod('NtClose').Invoke($null, ([object[]]@($directoryHandle)))
            $result.close_status = Status-Hex $closeStatus
            $result.directory_handle_closed = ($closeStatus -eq 0)
            if ($closeStatus -ne 0) { $result.gaps += 'NtClose did not report success' }
        } catch { $result.gaps += ('NtClose: ' + $_.Exception.Message) }
    }
    foreach ($allocation in @($queryBuffer, $attributes, $unicodeString, $nameBuffer)) {
        if ($allocation -ne [IntPtr]::Zero) { [Runtime.InteropServices.Marshal]::FreeHGlobal($allocation) }
    }
    $result.matches = $hits.ToArray()
    $result.elapsed_seconds = $watch.Elapsed.TotalSeconds
    $result.audit_success = ($result.query_completed -and $result.directory_handle_closed -and $result.gaps.Count -eq 0)
    $result | ConvertTo-Json -Depth 6 -Compress
    $process.Dispose()
}
if (-not $result.audit_success) { exit 2 }
