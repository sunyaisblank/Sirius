param([Parameter(Mandatory=$true)][string]$Config)
$ErrorActionPreference='Stop'
$settings=Get-Content -Raw -LiteralPath $Config | ConvertFrom-Json
if (Test-Path -LiteralPath $settings.guard_output) { throw 'Query guard output already exists' }
$folder=(New-Item -ItemType Directory -Path $settings.guard_output).FullName
$temporary=(New-Item -ItemType Directory -Path (Join-Path $folder 'temp')).FullName
$env:TEMP=$temporary;$env:TMP=$temporary
$self=Get-Process -Id $PID;$selfHandle=$self.Handle
$record=[ordered]@{schema='sirius-native-shader-query-owner-v1';owner_pid=$PID;owner_birth=$self.StartTime.ToUniversalTime().ToString('o');source_revision=$settings.source_revision;deadline_seconds=120;rss_guard_bytes=4294967296;scope='One live-driver compile/statistics query,zero dispatches;sampled childWorkingSet only,not GPU/aggregate memory or numerical/performance qualification';passed=$false;stop_reason=$null;errors=@()}
$record | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $folder 'owner.json') -Encoding UTF8
Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
public static class SiriusQueryExit {
 [DllImport("kernel32.dll", SetLastError=true)]
 [return: MarshalAs(UnmanagedType.Bool)]
 private static extern bool GetExitCodeProcess(IntPtr handle, out uint code);
 public static uint Read(IntPtr handle) {
  uint code;
  if(!GetExitCodeProcess(handle,out code))throw new System.ComponentModel.Win32Exception(Marshal.GetLastWin32Error());
  if(code==259)throw new Exception("Owned query process is still active");
  return code;
 }
}
'@
$process=$null;$samples=New-Object 'System.Collections.Generic.List[object]';$watch=[Diagnostics.Stopwatch]::StartNew()
try {
 $arguments=@('-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',('"'+$settings.query_script+'"'),'-Config',('"'+$Config+'"'))
 $process=Start-Process -FilePath (Join-Path $PSHOME 'powershell.exe') -ArgumentList $arguments -NoNewWindow -PassThru -RedirectStandardOutput (Join-Path $folder 'stdout.log') -RedirectStandardError (Join-Path $folder 'stderr.log')
 $handle=$process.Handle;$record.child_pid=$process.Id;$record.child_birth=$process.StartTime.ToUniversalTime().ToString('o')
 $record | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $folder 'owner.json') -Encoding UTF8
 while(!$process.HasExited) {
  $process.Refresh();if($process.HasExited){break}
  $samples.Add([pscustomobject]@{seconds=$watch.Elapsed.TotalSeconds;working_set_bytes=$process.WorkingSet64})
  if($process.WorkingSet64 -gt 4294967296){$record.stop_reason='sampled_rss'}
  if($watch.Elapsed.TotalSeconds -gt 120){$record.stop_reason='elapsed'}
  if($record.stop_reason){$process.Kill();break}
  Start-Sleep -Milliseconds 100
 }
 if(!$process.WaitForExit(5000)){throw 'Query process did not exit'}
 $record.returncode=[SiriusQueryExit]::Read($handle)
 if($record.returncode -ne 0 -or $record.stop_reason){throw 'Query child nonpass'}
 $record.passed=$true
} catch {
 $record.errors+=@($_.Exception.Message)
} finally {
 if($process -and !$process.HasExited){$process.Kill();if(!$process.WaitForExit(5000)){$record.errors+=@('Owned query child did not stop')}}
 $record.child_absent=(!$process -or $process.HasExited)
 $record.elapsed_seconds=$watch.Elapsed.TotalSeconds
 $record.peak_sampled_child_rss_bytes=if($samples.Count){($samples.ToArray()|Measure-Object -Property working_set_bytes -Maximum).Maximum}else{$null}
 $record.sample_count=$samples.Count
 ConvertTo-Json -InputObject $samples.ToArray() -Depth 5 | Set-Content -LiteralPath (Join-Path $folder 'samples.json') -Encoding UTF8
 $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $folder 'owner.json') -Encoding UTF8
}
if(!$record.passed -or !$record.child_absent -or $record.errors.Count){exit 1}
