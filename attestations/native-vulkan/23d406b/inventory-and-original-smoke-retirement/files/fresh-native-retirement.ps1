$ErrorActionPreference='Stop'
$records = [Text.Encoding]::UTF8.GetString([Convert]::FromBase64String('W3sia2luZCI6ICJpbnZlbnRvcnktb3duZXIiLCAicGlkIjogMjIzNTYsICJiaXJ0aCI6ICIyMDI2LTEwLTE1VDA4OjIyOjMwLjczMDgxMjZaIn0sIHsia2luZCI6ICJpbnZlbnRvcnktdGVzdCIsICJwaWQiOiAxODIxMiwgImJpcnRoIjogIjIwMjYtMTAtMTVUMDg6MjI6NDIuNjEwOTcwNVoifSwgeyJraW5kIjogImludmVudG9yeS10ZXJtaW5hbC1hdWRpdG9yIiwgInBpZCI6IDEwODgsICJiaXJ0aCI6ICIyMDI2LTEwLTE1VDA4OjIzOjIwLjMzMjk2MjdaIn0sIHsia2luZCI6ICJzbW9rZS1vd25lciIsICJwaWQiOiAxOTI3NiwgImJpcnRoIjogIjIwMjYtMTAtMTVUMDg6MjQ6MDkuODM5NjIzMVoifSwgeyJraW5kIjogInNtb2tlLXRlc3QiLCAicGlkIjogMjAyNjAsICJiaXJ0aCI6ICIyMDI2LTEwLTE1VDA4OjI0OjEyLjgyNDE1NzlaIn0sIHsia2luZCI6ICJzbW9rZS10ZXJtaW5hbC1hdWRpdG9yIiwgInBpZCI6IDE3MjgwLCAiYmlydGgiOiAiMjAyNi0xMC0xNVQwODoyNTowOS43NDM1NTMzWiJ9XQ==')) | ConvertFrom-Json
$rows = @()
foreach ($row in $records) {
 $p = Get-Process -Id ([int]$row.pid) -ErrorAction SilentlyContinue
 $same = $false
 $current = $null
 if ($p) { $current = $p.StartTime.ToUniversalTime().ToString('o'); $same = ($current -eq $row.birth) }
 $rows += [PSCustomObject]@{kind=$row.kind; pid=$row.pid; recorded_birth=$row.birth; current_birth=$current; recorded_birth_absent=(!$same)}
 if ($same) { throw 'Recorded owned native birth remains live' }
}
[PSCustomObject]@{schema='sirius-native-fresh-retirement-observation-v1'; observed_utc=[DateTime]::UtcNow.ToString('o'); passed=$true; recorded_births=$rows; scope='Read-only exact recorded native owner/child/auditor birth observations only; no provider launch, global process or hidden handle absence claim.'} | ConvertTo-Json -Depth 8
