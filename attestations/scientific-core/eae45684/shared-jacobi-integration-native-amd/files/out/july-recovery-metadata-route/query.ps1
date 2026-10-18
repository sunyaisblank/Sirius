$ErrorActionPreference='Stop'
[Console]::OutputEncoding=New-Object System.Text.UTF8Encoding($false)
$pattern='(?i)sirius|dngr|july|renderer|forensic|comparator|kerr|schwarzschild|alcubierre|wormhole|ellis|imax'
$result=[ordered]@{schema='sirius-july-recovery-metadata-route-v1';utc=[DateTime]::UtcNow.ToString('o');pid=$PID;process_start_utc=(Get-Process -Id $PID).StartTime.ToUniversalTime().ToString('o');user=[Security.Principal.WindowsIdentity]::GetCurrent().Name;pattern=$pattern;recent=$null;recycle=@();file_history=@();shadow_copies=@();gaps=@();read_only_metadata=$true}
try {
 $recent=[Environment]::GetFolderPath('Recent');$links=@(Get-ChildItem -LiteralPath $recent -Filter '*.lnk' -File -Force | Sort-Object LastWriteTimeUtc -Descending);$selected=@($links | Select-Object -First 1000);$shell=New-Object -ComObject WScript.Shell;$matches=@();$failures=0
 foreach($item in $selected){try{$link=$shell.CreateShortcut($item.FullName);if(($item.Name+' '+$link.TargetPath) -match $pattern){$matches+=@{link_name=$item.Name;link_mtime_utc=$item.LastWriteTimeUtc.ToString('o');target=$link.TargetPath;target_exists=[IO.File]::Exists($link.TargetPath) -or [IO.Directory]::Exists($link.TargetPath)}}}catch{$failures++}}
 $result.recent=@{path=$recent;total_links=$links.Count;inspected=$selected.Count;limit=1000;failures=$failures;matches=$matches}
} catch {$result.gaps+=('Recent: '+$_.Exception.Message)}
$sid=[Security.Principal.WindowsIdentity]::GetCurrent().User.Value
foreach($drive in [IO.DriveInfo]::GetDrives()){
 if($drive.DriveType -ne [IO.DriveType]::Fixed){continue}
 $base=Join-Path $drive.Name ('$Recycle.Bin\'+$sid)
 try{
  if(-not [IO.Directory]::Exists($base)){$result.recycle+=@{root=$base;exists=$false};continue}
  $items=@(Get-ChildItem -LiteralPath $base -Filter '$I*' -File -Force | Sort-Object LastWriteTimeUtc -Descending);$selected=@($items | Select-Object -First 2000);$matches=@();$failures=0
  foreach($item in $selected){try{
   if($item.Length -gt 65536 -or $item.Length -lt 24){$failures++;continue}
   $data=[IO.File]::ReadAllBytes($item.FullName);$version=[BitConverter]::ToInt64($data,0)
   if($version -eq 1){$offset=24;$length=$data.Length-24}elseif($version -eq 2 -and $data.Length -ge 28){$offset=28;$length=$data.Length-28}else{$failures++;continue}
   $original=[Text.Encoding]::Unicode.GetString($data,$offset,$length).TrimEnd([char]0)
   if($original -match $pattern){$payload=Join-Path $base ('$R'+$item.Name.Substring(2));$matches+=@{metadata=$item.FullName;original_path=$original;original_bytes=[BitConverter]::ToInt64($data,8);deleted_utc=[DateTime]::FromFileTimeUtc([BitConverter]::ToInt64($data,16)).ToString('o');payload=$payload;payload_exists=[IO.File]::Exists($payload) -or [IO.Directory]::Exists($payload)}}
  }catch{$failures++}}
  $result.recycle+=@{root=$base;exists=$true;metadata_total=$items.Count;inspected=$selected.Count;limit=2000;failures=$failures;matches=$matches}
 }catch{$result.gaps+=('Recycle '+$base+': '+$_.Exception.Message)}
}
foreach($name in @('Config1.xml','Config2.xml')){
 $path=Join-Path $env:LOCALAPPDATA ('Microsoft\Windows\FileHistory\Configuration\'+$name)
 try{
  if([IO.File]::Exists($path)){
   [xml]$xml=[IO.File]::ReadAllText($path);$refs=@($xml.SelectNodes('//*[contains(translate(local-name(),"ABCDEFGHIJKLMNOPQRSTUVWXYZ","abcdefghijklmnopqrstuvwxyz"),"target")]') | ForEach-Object {$_.InnerText})
   $result.file_history+=@{path=$path;exists=$true;target_references=$refs}
  }else{$result.file_history+=@{path=$path;exists=$false}}
 }catch{$result.gaps+=('FileHistory: '+$_.Exception.Message)}
}
try{$result.shadow_copies=@(Get-CimInstance -ClassName Win32_ShadowCopy -OperationTimeoutSec 10 | Select-Object ID,DeviceObject,InstallDate,VolumeName)}catch{$result.gaps+=('ShadowCopy: '+$_.Exception.Message)}
$result | ConvertTo-Json -Depth 8 -Compress

