$ErrorActionPreference='Stop'
$births=Get-Content -Raw -LiteralPath '\\wsl.localhost\Ubuntu\home\astra\.project\Sirius\out\native-kerr-context-preservation\native-births.json'|ConvertFrom-Json
$rows=@(foreach($b in $births){$p=Get-Process -Id ([int]$b.pid) -ErrorAction SilentlyContinue;$same=$false;if($p){$same=$p.StartTime.ToUniversalTime().ToString('o') -eq $b.start};if($same){throw 'Recorded native birth still live'};[ordered]@{case=$b.case;kind=$b.kind;pid=$b.pid;start=$b.start;recorded_birth_absent=(!$same)}})
[ordered]@{schema='sirius-recorded-native-birth-retirement-check-v1';rows=$rows;scope='Read-only exact PID/start check, no scientific launch/global absence claim'}|ConvertTo-Json -Depth 6
