"""Adapt existing bounded native ownership for a read-only feature query."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,subprocess
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def win(p):return subprocess.check_output(['wslpath','-w',str(p.resolve())],text=True).strip()
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def save(name,text):
 p=W/name;assert not p.exists();p.write_text(text)
owner=(W/'prior-owner.ps1').read_text()
owner=owner.replace("[int]$Seconds=320,[long]$RssBytes=4294967296", "[int]$Seconds=45,[long]$RssBytes=536870912")
anchor="Add-Type -TypeDefinition @'"
bootstrap=r"""$ownerProcess=Get-Process -Id $PID
$ownerHandle=$ownerProcess.Handle
$ownerStart=$ownerProcess.StartTime.ToUniversalTime().ToString('o')
$bootstrap=[ordered]@{owner_pid=$PID;owner_start_time_utc=$ownerStart;expected_test_executable='C:\Windows\System32\WindowsPowerShell\v1.0\powershell.exe';scope='Early exact native observer birth before helper compilation'}
$bootstrap | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $TemporaryDirectory 'bootstrap.json') -Encoding UTF8
"""
assert owner.count(anchor)==1;owner=owner.replace(anchor,bootstrap+anchor)
owner=owner.replace("$xml=Join-Path $Output 'gtest.xml'", "")
owner=owner.replace("$arguments=@('-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',(Join-Path $Stage 'child.ps1'))", "$arguments=@('-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',(Join-Path $Stage 'clock-child.ps1'),'-Stage',$Stage,'-OwnerTemporaryDirectory',$TemporaryDirectory)")
owner=owner.replace("$record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8\n while", "$record | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $Output 'owner.json') -Encoding UTF8\n [IO.File]::WriteAllText((Join-Path $Stage 'START.json'),(@{owner_pid=$PID;owner_start_time_utc=$ownerStart;child_pid=$process.Id;child_start_time_utc=$record.test_start_time_utc}|ConvertTo-Json),[Text.UTF8Encoding]::new($false))\n while")
owner=owner.replace("(Join-Path $Stage 'census.json')", "(Join-Path $Stage 'clock-features.json')")
owner=owner.replace("'Incomplete bounded extension census'", "'Incomplete bounded clock-feature query'")
owner=owner.replace("'3045067d4c29da744059798e175b97942a319f84'", "'6be3132b7023f9d6497ab8f67dc779cf8eb4a80c'")
owner=owner.replace("scope='Unadopted exact-module native extension census; independent numerical acceptance remains offline. Diagnostic RSS/time guard invalidates incomplete runs; no scientific/device budget changes or full runtime/release qualification claim.'", "scope='Read-only exact native shader-clock feature query; no logical device, pipeline, dispatch or scientific/release qualification. Fixed 45s/512MiB diagnostic guard.'")
owner=owner.replace("$busy=@(Get-Process | Where-Object {$_.ProcessName -match '^sirius'})\nif ($busy.Count) { throw 'A native Sirius consumer is already live' }\n", "")
save('clock-owner.ps1',owner)
child=r'''param([Parameter(Mandatory=$true)][string]$Stage,[Parameter(Mandatory=$true)][string]$OwnerTemporaryDirectory)
$ErrorActionPreference='Stop'
[Console]::OutputEncoding=[Text.UTF8Encoding]::new($false)
$self=Get-Process -Id $PID;$held=$self.Handle;$birth=$self.StartTime.ToUniversalTime().ToString('o')
try {
 $nativeTemp=New-Item -ItemType Directory -Path (Join-Path $Stage 'native-temp')
 $env:TEMP=$nativeTemp.FullName;$env:TMP=$env:TEMP
 $bootstrap=[ordered]@{child_pid=$PID;child_start_time_utc=$birth;scope='Exact child birth held before compiler/provider execution'}
 [IO.File]::WriteAllText((Join-Path $Stage 'child-bootstrap.json'),($bootstrap|ConvertTo-Json),[Text.UTF8Encoding]::new($false))
 $owner=Get-Content -Raw -LiteralPath (Join-Path $OwnerTemporaryDirectory 'bootstrap.json')|ConvertFrom-Json
 $clock=[Diagnostics.Stopwatch]::StartNew()
 while(!(Test-Path -LiteralPath (Join-Path $Stage 'START.json'))) {
  if($clock.Elapsed.TotalSeconds -gt 20){throw 'Parent START barrier exceeded20s'}
  $parent=Get-Process -Id $owner.owner_pid -ErrorAction SilentlyContinue
  if(!$parent -or $parent.StartTime.ToUniversalTime().ToString('o') -ne $owner.owner_start_time_utc){throw 'Recorded parent absent before START'}
  Start-Sleep -Milliseconds 20
 }
 $start=Get-Content -Raw -LiteralPath (Join-Path $Stage 'START.json')|ConvertFrom-Json
 if($start.child_pid -ne $PID -or $start.child_start_time_utc -ne $birth -or $start.owner_pid -ne $owner.owner_pid -or $start.owner_start_time_utc -ne $owner.owner_start_time_utc){throw 'START does not bind exact owner/child births'}
 Add-Type -TypeDefinition ([IO.File]::ReadAllText((Join-Path $Stage 'native_clock_features.cs')))
 $result=[SiriusNativeClockFeatures]::Query('C:\Windows\System32\vulkan-1.dll',[int]$owner.owner_pid,[string]$owner.owner_start_time_utc)
 [IO.File]::WriteAllText((Join-Path $Stage 'clock-features.json'),($result|ConvertTo-Json -Depth 10),[Text.UTF8Encoding]::new($false))
 [IO.File]::WriteAllText((Join-Path $Stage 'run/driver-modules.json'),(ConvertTo-Json -InputObject @($result['loaded_modules']) -Depth 8),[Text.UTF8Encoding]::new($false))
 [Console]::WriteLine('Native clock feature query completed')
 [Environment]::Exit(0)
}catch{[Console]::Error.WriteLine($_.Exception.ToString());[Environment]::Exit(1)}
'''
save('clock-child.ps1',child)
inputs={str(p.relative_to(ROOT)):seal(p) for p in (Path(__file__),W/'prior-owner.ps1',W/'native_clock_features.cs')}
save('owner-preparation.json',json.dumps({'input_seals':inputs,'owner':seal(W/'clock-owner.ps1'),'child':seal(W/'clock-child.ps1'),'windows_stage':win(W),'scope':'Source preparation only; early birth, held child barrier, fixed guarded read-only query and task-local compiler temporary directories.'},indent=2)+'\n')
print(json.dumps({'owner':seal(W/'clock-owner.ps1'),'child':seal(W/'clock-child.ps1')}))
