param([Parameter(Mandatory=$true)][string]$Stage,[Parameter(Mandatory=$true)][string]$NativeStage,[Parameter(Mandatory=$true)][string]$Receipt,[Parameter(Mandatory=$true)][string]$Bootstrap)
$ErrorActionPreference='Stop'
$self=Get-Process -Id $PID;$held=$self.Handle
@{owner_pid=$PID;owner_start_time_utc=$self.StartTime.ToUniversalTime().ToString('o')}|ConvertTo-Json|Set-Content -LiteralPath $Bootstrap -Encoding UTF8
$record=[ordered]@{kind='finite-native-clock-control-retirement';observer_pid=$PID;observer_birth=$self.StartTime.ToUniversalTime().ToString('o');targets=@();matches=@();errors=@();pass=$false;stage=$Stage;native_stage=$NativeStage;scope='Exact recorded owner/child births and accessible command-line references to this WORK and native stage. No hidden/global/scheduled-consumer absence claim.'}
try {
 $targets=@()
 $actions=(Get-Content -Raw -LiteralPath (Join-Path $Stage 'native-cases.json')|ConvertFrom-Json).actions
 if($actions.Count -ne 16){throw 'Expected all16trialactions'}
 foreach($action in $actions){
  $owner=Get-Content -Raw -LiteralPath (Join-Path $Stage ($action+'/run/owner.json'))|ConvertFrom-Json
  $targets+=@(@{kind=$action+' owner';pid=[int]$owner.owner_pid;birth=$owner.owner_start_time_utc},@{kind=$action+' child';pid=[int]$owner.test_pid;birth=$owner.test_start_time_utc})
 }
 foreach($target in $targets){
  $p=Get-Process -Id $target.pid -ErrorAction SilentlyContinue
  if($p){$handle=$p.Handle}
  $target.recorded_birth_absent=(!$p -or $p.StartTime.ToUniversalTime().ToString('o') -ne $target.birth)
  $record.targets+=@($target)
  if(!$target.recorded_birth_absent){throw 'Recorded native control birth remains live'}
 }
 $processes=@(Get-CimInstance Win32_Process -OperationTimeoutSec 20)
 $record.visible_processes=$processes.Count
 $record.hidden_field_records=@($processes|Where-Object {!$_.CommandLine -or !$_.ExecutablePath}).Count
 $record.matches=@($processes|Where-Object {[int]$_.ProcessId -ne $PID -and $_.CommandLine -and ($_.CommandLine.IndexOf($Stage,[StringComparison]::OrdinalIgnoreCase) -ge 0 -or $_.CommandLine.IndexOf($NativeStage,[StringComparison]::OrdinalIgnoreCase) -ge 0)}|Select-Object ProcessId,ParentProcessId,CreationDate,Name,ExecutablePath,CommandLine)
 if($record.matches.Count){throw 'Accessible consumer still refers to task WORK'}
 $record.pass=$true
}catch{$record.errors+=@($_.Exception.Message)}
$record|ConvertTo-Json -Depth 8|Set-Content -LiteralPath $Receipt -Encoding UTF8
if(!$record.pass){exit 1}
