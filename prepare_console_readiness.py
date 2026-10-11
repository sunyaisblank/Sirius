"""Derive a tools-only real-console readiness diagnostic from reviewed owners."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,ast
w=Path(__file__).resolve().parent
def derive(name,expected):
 p=w/name;b=p.read_bytes();assert hashlib.sha256(b).hexdigest()==expected;return b.decode()
ps=derive('service_info_loopback.ps1','d3cdc894f171061e58fe5a73140803d555ebb90a17bd2dd47bb010b7b7785b7a')
def change(old,new):
 global ps
 assert ps.count(old)==1,old;ps=ps.replace(old,new)
ps=ps.replace('service-info-loopback','console-readiness')
change("$argument=@('--system-info','--verbose','--process','sirius_backend_tests.exe','--remote-host','127.0.0.1','--remote-port','27300')", """$argument=@('--mode=profiling','--verbose','--process=sirius_backend_tests.exe','--remote-host=127.0.0.1','--remote-port=27300','--rgp-auto-capture=dispatch:43:2','--rgp-capture-mode=dispatch','--rgp-render-op-count=2','--rgp-instruction-tracing','--rgp-counter-collection','--rgp-sqtt-buffer-size=default',('"--output='+(Join-Path $Output 'unused.rgp')+'"'))""")
change("Read-only system information only; no scientific workload or capture.","Owned real-console readiness inspection only; no scientific application or actual trace.")
change("'@\n$existing=", "'@\nAdd-Type -TypeDefinition (Get-Content -Raw -LiteralPath (Join-Path $PSScriptRoot 'owned_console.cs'))\n$existing=")
change("Pinned owned service and system-info client only. Deliberate owned service stop after information request; no application capture, clock/driver experiment or scientific qualification.","Unchanged auto-capture arguments in an owned real console. Require actual enabled/waiting screen text, then deliberately stop tools; no application/dispatch/actual trace or qualification.")
change(" $process=Start-Process -FilePath $exe -ArgumentList $argument -WorkingDirectory $Tools -NoNewWindow -RedirectStandardOutput (Join-Path $Output 'stdout.log') -RedirectStandardError (Join-Path $Output 'stderr.log') -PassThru", """ $consolePid=[SiriusOwnedConsole]::StartSuspended($exe,($argument -join ' '),$Tools)
 $process=Get-Process -Id $consolePid""")
change(" $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8\n while (!$process.HasExited)", """ $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 $attached=$false;$record.console_ready_observed=$false
 $attach=[Diagnostics.Stopwatch]::StartNew()
 while (!$attached) {
  $service.Refresh();$process.Refresh()
  $peak=[Math]::Max($peak,[Math]::Max($process.WorkingSet64,$service.WorkingSet64))
  if ($service.HasExited -or $process.HasExited -or $attach.Elapsed.TotalSeconds -gt 3 -or $peak -gt 1073741824) {throw 'Owned suspended-console attachment/lifetime/RSS bound exceeded'}
  $attached=[SiriusOwnedConsole]::TryAttach([uint32]$process.Id)
  if (!$attached) {Start-Sleep -Milliseconds 50}
 }
 $initial=[SiriusOwnedConsole]::Read([uint32]$PID,[uint32]$process.Id)
 $initial | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'initial-screen.json') -Encoding UTF8
 $record.console_geometry_initialized_before_resume=$true
 $record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8
 [SiriusOwnedConsole]::Resume()
 while (!$process.HasExited)""")
change("  Start-Sleep -Milliseconds 50", """  if (!$attached) {$attached=[SiriusOwnedConsole]::TryAttach([uint32]$process.Id)}
  if ($attached) {
   $snapshot=[SiriusOwnedConsole]::Read([uint32]$PID,[uint32]$process.Id)
   if ($snapshot.text.Contains('Profiling (RGP) feature enabled (auto-capture enabled)') -and $snapshot.text.Contains('Waiting for GPU application to connect...')) {
    $snapshot | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'ready-screen.json') -Encoding UTF8
    $snapshot.text | Set-Content -LiteralPath (Join-Path $Output 'stdout.log') -Encoding UTF8
    $record.console_ready_observed=$true;$record.readiness_elapsed_seconds=$watch.Elapsed.TotalSeconds
    $record.observation_only_intentional_tool_stop=$true
    $process.Kill();if (!$process.WaitForExit(5000)) {throw 'Observed owned console tool did not stop'}
    break
   }
  }
  Start-Sleep -Milliseconds 50""")
change(" $record.exit_code_source='GetExitCodeProcess on preacquired owned handle after WaitForExit'", """ $record.exit_code_source='GetExitCodeProcess on preacquired owned handle after WaitForExit'
 if (!$record.console_ready_observed -or !$record.observation_only_intentional_tool_stop) {throw 'Console readiness was not observed; no scientific launch'}""")
change(" $record.owned_service_absent=($null -eq $service)", """ $record.console_detached=$false
 try {[SiriusOwnedConsole]::Detach();$record.console_detached=$true} catch {$errors.Add('console-detach: '+$_.Exception.Message)}
 $record.owned_service_absent=($null -eq $service)""")
change("if (!$record.owned_service_absent -or", "if (!$record.console_detached -or !$record.owned_service_absent -or")
change("if ($record.status -ne 'completed' -or $record.returncode -ne 0)", "if ($record.status -ne 'completed' -or !$record.console_ready_observed -or !$record.observation_only_intentional_tool_stop -or !$record.console_detached)")
with (w/'console_readiness.ps1').open('x') as f:f.write(ps)
py=derive('run_service_info_loopback.py','212d96484073badd439837fb940307f5ec28619166c5c9242d50c94b5fd648e0')
py=py.replace('service-info-loopback','console-readiness').replace('service_info_loopback.ps1','console_readiness.ps1')
py=py.replace("paths += [ROOT/v['path'] for v in extraction['files']]", "paths += [WORK/'owned_console.cs',WORK/'prepare_console_readiness.py']\npaths += [ROOT/v['path'] for v in extraction['files']]")
py=py.replace("o['status']=='completed' and o['returncode']==0", "o['status']=='completed' and o['console_ready_observed'] and o['observation_only_intentional_tool_stop'] and o['console_detached'] and o['console_geometry_initialized_before_resume']")
needle="record['stdout_sha256']=sha(output/'stdout.log');record['stderr_sha256']=sha(output/'stderr.log')"
assert py.count(needle)==1
py=py.replace(needle, "record['stdout_sha256']=sha(output/'stdout.log')\n  screen=read(output/'ready-screen.json');assert screen['width']==1024 and screen['height']==2048 and screen['cursor_y']<2047 and screen['rows_read']==screen['cursor_y']+1 and len(screen['cells'])==screen['width']*screen['rows_read']\n  assert set(screen['attached_pids'])=={o['owner_pid'],o['test_pid']}\n  assert 'Profiling (RGP) feature enabled (auto-capture enabled)' in screen['text'] and 'Waiting for GPU application to connect...' in screen['text']\n  assert not (output/'unused.rgp').exists()\n  record['ready_screen_sha256']=sha(output/'ready-screen.json')")
py=py.replace('Owned RDS and read-only system-info request only. No scientific application or capture;', 'Owned RDS and unchanged auto-capture client in real console; actual readiness observation followed by deliberate stop. No scientific application, dispatch or actual trace;')
ast.parse(py)
with (w/'run_console_readiness.py').open('x') as f:f.write(py)
for n in ('console_readiness.ps1','run_console_readiness.py','owned_console.cs',Path(__file__).name):
 p=w/n;print(n,p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest())
