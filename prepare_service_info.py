"""Derive a two-process, read-only RDS connectivity diagnostic."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib
work = Path(__file__).resolve().parent
ps = (work / 'tool_help.ps1').read_text(); bridge = (work / 'run_tool_help.py').read_text()
assert hashlib.sha256(ps.encode()).hexdigest() == 'b5f166858f531c7a0edf565ad7cd22385078f5f6f3c95f75501e2f64c3bdc5cd'
assert hashlib.sha256(bridge.encode()).hexdigest() == 'a63df2d251fa4c7f16864daa7d9f43e3c498807acc9ca9dfa75292d35b521c05'
def change(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new)
ps = change(ps, "[ValidateSet('panel-help','panel-version','profiler-help')]", "[ValidateSet('service-info')]")
ps = change(ps, '[string]$ExpectedHash,', '[string]$ExpectedHash,[Parameter(Mandatory=$true)][string]$ExpectedServiceHash,')
ps = change(ps, "$tool=if ($Action -eq 'profiler-help') {'RadeonGPUProfiler.exe'} else {'RadeonDeveloperPanelCLI.exe'}\n$argument=if ($Action -eq 'panel-version') {'--version'} else {'--help'}", "$tool='RadeonDeveloperPanelCLI.exe'\n$argument=@('--system-info','--verbose','--process','sirius_backend_tests.exe')")
ps = change(ps, '$exe=Join-Path $Tools $tool', "$exe=Join-Path $Tools $tool\n$serviceExe=Join-Path $Tools 'RadeonDeveloperServiceCLI.exe'")
ps = change(ps, 'expected_test_executable=$exe;', 'expected_test_executable=$exe;owned_child_executables=@($exe,$serviceExe);')
ps = change(ps, "if ($existing.Count) { throw 'Existing developer/profiler process: tool help declined' }", "if ($existing.Count -or @(Get-Process | Where-Object {$_.ProcessName -match '^sirius'}).Count) { throw 'Existing developer/profiler/Sirius consumer: service diagnostic declined' }")
ps = change(ps, "$record=[ordered]", """$serviceHash=(Get-FileHash -LiteralPath $serviceExe -Algorithm SHA256).Hash.ToLowerInvariant()
if ($serviceHash -ne $ExpectedServiceHash) {throw 'Pinned service mismatch'}
$serviceSignature=Get-AuthenticodeSignature -LiteralPath $serviceExe
@{status=$serviceSignature.Status.ToString();subject=$serviceSignature.SignerCertificate.Subject;thumbprint=$serviceSignature.SignerCertificate.Thumbprint;file_version=(Get-Item -LiteralPath $serviceExe).VersionInfo.FileVersion} | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'service-signature.json') -Encoding UTF8
if ($serviceSignature.Status -ne 'Valid') {throw 'Official service signature did not verify'}
$record=[ordered]""")
ps = change(ps, "status='running';scope=", "status='running';scoped_processes=@();owned_service_absent=$false;scope=")
ps = change(ps, '$process=$null;$watch=', '$process=$null;$service=$null;$watch=')
ps = change(ps, "try {\n $process=Start-Process", """try {
 $service=Start-Process -FilePath $serviceExe -WorkingDirectory $Tools -NoNewWindow -RedirectStandardOutput (Join-Path $Output 'service-stdout.log') -RedirectStandardError (Join-Path $Output 'service-stderr.log') -PassThru
 $serviceHandle=$service.Handle
 $record.scoped_processes=@(@{kind='owned-RDS';pid=$service.Id;birth=$service.StartTime.ToUniversalTime().ToString('o');executable=$serviceExe;sha256=$serviceHash})
 $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 Start-Sleep -Milliseconds 1000
 if ($service.HasExited) {throw 'Owned service exited before the information client'}
 $process=Start-Process""")
ps = change(ps, '  $peak=[Math]::Max($peak,$process.WorkingSet64)', "  $service.Refresh();if ($service.HasExited) {throw 'Owned service exited during information request'}\n  $peak=[Math]::Max($peak,[Math]::Max($process.WorkingSet64,$service.WorkingSet64))")
ps = change(ps, '-or $process.WorkingSet64 -gt 1073741824)', '-or $process.WorkingSet64 -gt 1073741824 -or $service.WorkingSet64 -gt 1073741824)')
ps = change(ps, ' $record.elapsed_seconds=', """ $record.owned_service_absent=($null -eq $service)
 if ($service) {
  try {
   if (!$service.HasExited) {$record.intentional_service_stop=$true;$service.Kill();if (!$service.WaitForExit(5000)) {throw 'Owned service did not terminate'}}
   $record.service_returncode=[SiriusOwnedExit]::Read($serviceHandle)
   $record.owned_service_absent=$service.HasExited
  } catch {$errors.Add('service-stop: '+$_.Exception.Message)}
 }
 $record.service_unchanged=$false
 try {$record.service_unchanged=((Get-FileHash -LiteralPath $serviceExe -Algorithm SHA256).Hash.ToLowerInvariant() -eq $serviceHash)} catch {$errors.Add('service-posthash: '+$_.Exception.Message)}
 $record.elapsed_seconds=""")
ps = change(ps, 'if ($errors.Count -or !$record.owned_process_absent', 'if (!$record.owned_service_absent -or !$record.service_unchanged -or $errors.Count -or !$record.owned_process_absent')
ps = ps.replace('Early exact tool-help observer birth; no service/capture/scientific launch authorized.', 'Early exact service/client owner birth. Read-only system information only; no scientific workload or capture.')
ps = ps.replace('Pinned official executable help/version only; no service/capture/scientific result or developer-mode/clock/driver experiment authorized.', 'Pinned owned service and system-info client only. Deliberate owned service stop after information request; no application capture, clock/driver experiment or scientific qualification.')
bridge = change(bridge, "('panel-help','panel-version','profiler-help')", "('service-info',)")
bridge = change(bridge, "tool='RadeonGPUProfiler.exe' if action=='profiler-help' else 'RadeonDeveloperPanelCLI.exe'", "tool='RadeonDeveloperPanelCLI.exe'")
bridge = bridge.replace('tool-help-source-review.json','service-info-source-review.json').replace('tool_help.ps1','service_info.ps1').replace('emergency_native_cleanup.ps1','emergency_profile_cleanup.ps1')
bridge = change(bridge, "exe=matches[0];tools=exe.parent;output=", "exe=matches[0];tools=exe.parent;service=tools/'RadeonDeveloperServiceCLI.exe';output=")
bridge = change(bridge, "'-ExpectedHash',sha(exe),'-LiveRevision'", "'-ExpectedHash',sha(exe),'-ExpectedServiceHash',sha(service),'-LiveRevision'")
bridge = change(bridge, "and o['developer_census_completed']", "and o['owned_service_absent'] and o['service_unchanged'] and o['developer_census_completed']")
bridge = change(bridge, "assert read(output/'signature.json')['status']=='Valid'", "assert read(output/'signature.json')['status']=='Valid' and read(output/'service-signature.json')['status']=='Valid'")
bridge = bridge.replace('Tool help/version only. No service/science/device capture;', 'Owned RDS and read-only system-info request only. No scientific application or capture;')
emergency = (work/'emergency_native_cleanup.ps1').read_text()
emergency = change(emergency, "if ($owner.test_pid -and $owner.test_start_time_utc)", """foreach ($item in @($owner.scoped_processes)) {
   if (!$item.pid -or !$item.birth -or !$item.executable -or !(@($bootstrap.owned_child_executables) -contains $item.executable)) {throw 'Unbound scoped process receipt'}
   $targets.Add([ordered]@{kind=$item.kind;pid=[int]$item.pid;birth=$item.birth})
  }
  if ($owner.test_pid -and $owner.test_start_time_utc)""")
emergency = change(emergency, "$isTest=($item.ExecutablePath -and [String]::Equals($item.ExecutablePath,$bootstrap.expected_test_executable,[StringComparison]::OrdinalIgnoreCase))", "$isTest=($item.ExecutablePath -and @($bootstrap.owned_child_executables | Where-Object {[String]::Equals($_,$item.ExecutablePath,[StringComparison]::OrdinalIgnoreCase)}).Count -gt 0)")
for name,data in [('service_info.ps1',ps),('run_service_info.py',bridge),('emergency_profile_cleanup.ps1',emergency)]:
    with (work/name).open('x') as f: f.write(data)
    print(name,len(data.encode()),hashlib.sha256(data.encode()).hexdigest())
