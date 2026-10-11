"""Derive an owned-console capture without changing the original scientific case."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
WORK=Path(__file__).resolve().parent
def replace(s,old,new):
    assert s.count(old)==1,(old[:100],s.count(old))
    return s.replace(old,new)
def save(name,s):
    p=WORK/name
    if p.exists():assert p.read_text()==s
    else:p.write_text(s)
cs=(WORK/'owned_console_requested.cs').read_text()
cs=replace(cs,'[Out]StringBuilder text,uint length','[Out,MarshalAs(UnmanagedType.LPArray,ArraySubType=UnmanagedType.U2,SizeParamIndex=2)]char[] text,uint length')
cs=replace(cs,'public int width,height,cursor_x,cursor_y,rows_read;','public int width,height,cursor_x,cursor_y,rows_read;\n  public uint requested_count,read_count;')
cs=replace(cs,'public static Snapshot Read(uint owner,uint child)','public static Snapshot Read(uint owner,uint child,uint science)')
cs=replace(cs,'else if(ids[i]!=child)','else if(ids[i]!=child && (science==0 || ids[i]!=science))')
cs=replace(cs,'var cells=new StringBuilder((int)length);uint read;','var cells=new char[(int)length];uint read;')
cs=replace(cs,'string raw=cells.ToString();var text=new StringBuilder();','string raw=new string(cells,0,(int)read);\n  BufferInfo after;if(!GetConsoleScreenBufferInfo(screen,out after))throw Error();\n  if(after.size.x!=info.size.x || after.size.y!=info.size.y || after.cursor.x!=info.cursor.x || after.cursor.y!=info.cursor.y)throw new Exception("Console changed during counted snapshot");\n  var text=new StringBuilder();')
cs=replace(cs,'rows_read=info.cursor.y+1,attached_pids=ids','rows_read=info.cursor.y+1,requested_count=length,read_count=read,attached_pids=ids')
save('owned_console_capture.cs',cs)
ps=(WORK/'run_capture_control.ps1').read_text()
ps=replace(ps,'[Parameter(Mandatory=$true)][string]$ServiceHash)','[Parameter(Mandatory=$true)][string]$ServiceHash,[Parameter(Mandatory=$true)][int]$RequiredConsoleWidth)')
ps=replace(ps,"'@\nforeach", "'@\nAdd-Type -TypeDefinition (Get-Content -Raw -LiteralPath (Join-Path $PSScriptRoot 'owned_console_capture.cs'))\nforeach")
ps=replace(ps,"$panelStdout=Join-Path $Output 'panel-stdout.log';$panelStderr=Join-Path $Output 'panel-stderr.log'","$panelScreenText=Join-Path $Output 'panel-screen.txt'")
ps=replace(ps,"$toolPeak=0\ntry", "$toolPeak=0;$consoleAttached=$false;$lastScreen=$null\ntry")
old=""" $panel=Start-Process -FilePath $panelExe -ArgumentList $panelArguments -WorkingDirectory $Tools -NoNewWindow -RedirectStandardOutput $panelStdout -RedirectStandardError $panelStderr -PassThru
 $panelHandle=$panel.Handle"""
new=""" [SiriusOwnedConsole]::ConfigureWidth($RequiredConsoleWidth)
 $record.required_console_width=$RequiredConsoleWidth
 $record.capture_output_kind='Counted owned-console screen cells; combined channels, not stdout byte history'
 $consolePid=[SiriusOwnedConsole]::StartSuspended($panelExe,($panelArguments -join ' '),$Tools)
 $panel=Get-Process -Id $consolePid
 $panelHandle=$panel.Handle"""
ps=replace(ps,old,new)
old=""" $ready=[Diagnostics.Stopwatch]::StartNew();$readyObserved=$false
 while (!$readyObserved)"""
new=""" [SiriusOwnedConsole]::Resume()
 $attach=[Diagnostics.Stopwatch]::StartNew()
 while (!$consoleAttached) {
  $service.Refresh();$panel.Refresh()
  $toolPeak=[Math]::Max($toolPeak,[Math]::Max($service.WorkingSet64,$panel.WorkingSet64))
  if ($service.HasExited -or $panel.HasExited -or $attach.Elapsed.TotalSeconds -gt 3 -or $toolPeak -gt 1073741824) {throw 'Owned running-console attachment/lifetime/RSS bound exceeded; no science launch'}
  $consoleAttached=[SiriusOwnedConsole]::TryAttach([uint32]$panel.Id)
  $record.last_console_attach_error=[SiriusOwnedConsole]::LastAttachError
  if (!$consoleAttached) {Start-Sleep -Milliseconds 50}
 }
 $lastScreen=[SiriusOwnedConsole]::Read([uint32]$PID,[uint32]$panel.Id,0)
 $lastScreen | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'initial-screen.json') -Encoding UTF8
 $record.console_geometry_verified_after_attachment=$true
 $ready=[Diagnostics.Stopwatch]::StartNew();$readyObserved=$false
 while (!$readyObserved)"""
ps=replace(ps,old,new)
ps=replace(ps,"  $text=Get-Content -Raw -LiteralPath $panelStdout\n  $readyObserved=($text -and $text.Contains('Waiting for GPU application to connect...'))", "  $lastScreen=[SiriusOwnedConsole]::Read([uint32]$PID,[uint32]$panel.Id,0)\n  $readyObserved=($lastScreen.text.Contains('Profiling (RGP) feature enabled (auto-capture enabled)') -and $lastScreen.text.Contains('Waiting for GPU application to connect...'))")
ps=replace(ps,' $record.capture_ready_observed_before_application=$true'," $lastScreen | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'ready-screen.json') -Encoding UTF8\n $record.capture_ready_observed_before_application=$true")
# The original application still owns its unchanged redirected stdout/XML.
# Its only allowed console membership is the exact held, recorded child PID.
ps=replace(ps," $record.capture_client_returncode=[SiriusOwnedExit]::Read($panelHandle)", " $lastScreen=[SiriusOwnedConsole]::Read([uint32]$PID,[uint32]$panel.Id,[uint32]$process.Id)\n $lastScreen | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'final-screen.json') -Encoding UTF8\n $lastScreen.text | Set-Content -LiteralPath $panelScreenText -Encoding UTF8\n $record.capture_client_returncode=[SiriusOwnedExit]::Read($panelHandle)")
ps=replace(ps," $record.owned_service_absent=($null -eq $service)\n if ($service)", " $record.console_detached=$false\n try {if ($lastScreen) {$lastScreen | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'last-screen.json') -Encoding UTF8}} catch {$cleanupErrors.Add('last-screen-receipt: '+$_.Exception.Message)}\n try {[SiriusOwnedConsole]::Detach();$record.console_detached=$true} catch {$cleanupErrors.Add('console-detach: '+$_.Exception.Message)}\n $record.owned_service_absent=($null -eq $service)\n if ($service)")
ps=replace(ps,' if (!$record.owned_capture_client_absent or'.replace(' or',' -or'),' if (!$record.console_detached -or !$record.owned_capture_client_absent -or')
ps=replace(ps,"-or !$record.capture_trace_sealed -or", "-or !$record.console_detached -or !$record.capture_trace_sealed -or")
save('run_console_capture_control.ps1',ps)
reader=(WORK/'verify_capture_control.py').read_text()
reader=replace(reader,"'developer_census_completed')","'developer_census_completed','console_detached','console_geometry_verified_after_attachment')")
reader=replace(reader,"/native-capture/original-modern-te.rgp'","/native-capture-console/original-modern-te.rgp'")
reader=replace(reader,"log=output/'panel-stdout.log'", "log=output/'panel-screen.txt'")
reader=replace(reader,"    text=log.read_text(encoding='utf-8-sig')", """    screens={}
    for name in ('initial','ready','final'):
        s=read(output/(name+'-screen.json'));screens[name]=s
        assert o['required_console_width']<=s['width']<=4096 and 2<=s['height']<=32767
        assert 0<=s['cursor_x']<s['width'] and 0<=s['cursor_y']<s['height']-1 and s['rows_read']==s['cursor_y']+1
        count=s['width']*s['rows_read'];assert 0<count<=2097152
        assert s['requested_count']==s['read_count']==count==len(s['cells'])
        lines=[s['cells'][i:i+s['width']] for i in range(0,count,s['width'])]
        assert all(line[-1]==' ' for line in lines)
        assert s['text']==''.join(line.rstrip(' ')+'\\r\\n' for line in lines)
        assert o['owner_pid'] in s['attached_pids'] and len(set(s['attached_pids']))==len(s['attached_pids'])
        allowed={o['owner_pid'],next(v['pid'] for v in o['scoped_processes'] if v['kind']=='owned-RDP-capture')}
        if name=='final':allowed.add(o['test_pid'])
        assert set(s['attached_pids'])<=allowed
        if name!='final':assert set(s['attached_pids'])==allowed
    assert len({(s['width'],s['height']) for s in screens.values()})==1
    assert screens['initial']['cursor_y']<=screens['ready']['cursor_y']<=screens['final']['cursor_y']
    assert 'Profiling (RGP) feature enabled (auto-capture enabled)' in screens['ready']['text']
    assert 'Waiting for GPU application to connect...' in screens['ready']['text']
    text=log.read_bytes().decode('utf-8-sig')
    assert text==screens['final']['text']+'\\r\\n'
    text=screens['final']['text']""")
save('verify_console_capture_control.py',reader)
bridge=(WORK/'run_native_capture.py').read_text()
bridge=replace(bridge,"('native-capture',)","('native-capture-console',)")
assert bridge.count("WORK/'native-capture-source-review.json'")==2
bridge=bridge.replace("WORK/'native-capture-source-review.json'", "WORK/'native-capture-console-source-review.json'")
bridge=bridge.replace("WORK/'run_capture_control.ps1'", "WORK/'run_console_capture_control.ps1'")
bridge=bridge.replace("WORK/'verify_capture_control.py'", "WORK/'verify_console_capture_control.py'")
bridge=replace(bridge,'from verify_capture_control import','from verify_console_capture_control import')
bridge=replace(bridge,"paths += [WORK/'prepare_capture_controller.py'", "paths += [WORK/'owned_console_capture.cs',WORK/'prepare_console_capture.py',WORK/'verify_requested_prefix.py',WORK/'console-requested-prefix-result.json',WORK/'prepare_capture_controller.py'")
bridge=replace(bridge,"ps='/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'", """science=win(exe);future_trace=win(output/'original-modern-te.rgp')
required_width=max(len('Connecting to process: '+science+' (PID: '+str(2**32-1)+')'),len('Trace saved to: '+future_trace),len('Profiling (RGP) feature enabled (auto-capture enabled)'),len('Waiting for GPU application to connect...'))+1
assert read(WORK/'console-requested-prefix-result.json')['bounded_readiness_prefix_verified']
ps='/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'""")
bridge=replace(bridge,"'-ServiceHash',sha(service)]", "'-ServiceHash',sha(service),'-RequiredConsoleWidth',str(required_width)]")
save('run_native_console_capture.py',bridge)
