$ErrorActionPreference='Stop'
$p=Get-Process -Id $PID;$h=$p.Handle
$base=Split-Path -Parent $MyInvocation.MyCommand.Path
@{owner_pid=$PID;owner_start_time_utc=$p.StartTime.ToUniversalTime().ToString('o')}|ConvertTo-Json|Set-Content -LiteralPath (Join-Path $base 'search-index-bootstrap.json') -Encoding UTF8
$query="SELECT TOP 60 System.ItemPathDisplay,System.ItemNameDisplay,System.DateModified,System.Size FROM SystemIndex WHERE (SCOPE='file:C:/Users/alarisadmin/Desktop' OR SCOPE='file:C:/Users/alarisadmin/Documents' OR SCOPE='file:C:/Users/alarisadmin/Downloads' OR SCOPE='file:C:/Users/alarisadmin/OneDrive') AND (System.FileName LIKE '%Sirius%' OR System.FileName LIKE '%render%' OR System.FileName LIKE '%comparator%' OR System.FileName LIKE '%tape%')"
$connection=$null;$command=$null;$reader=$null
try {
 $connection=New-Object System.Data.OleDb.OleDbConnection("Provider=Search.CollatorDSO;Extended Properties='Application=Windows';")
 $connection.Open();$command=$connection.CreateCommand();$command.CommandText=$query;$command.CommandTimeout=15
 $reader=$command.ExecuteReader();$indexRows=New-Object 'System.Collections.Generic.List[object]'
 while($reader.Read()){
  $row=[ordered]@{}
  for($i=0;$i -lt $reader.FieldCount;$i++){ $value=$reader.GetValue($i);if($value -is [DBNull]){$value=$null};if($value -is [DateTime]){$value=$value.ToUniversalTime().ToString('o')};$row[$reader.GetName($i)]=$value }
  $indexRows.Add($row)
 }
 @{utc=[DateTime]::UtcNow.ToString('o');query=$query;rows=$indexRows.ToArray();truncated_at_limit=($indexRows.Count -eq 60);scope='Read-only Windows Search filename metadata in four current-user roots; no file opening/hydration/restore/content scan. Index omissions/staleness and unindexed storage remain outside this route.'}|ConvertTo-Json -Depth 8
} catch {
 @{utc=[DateTime]::UtcNow.ToString('o');query=$query;error=$_.Exception.Message;scope='Windows Search metadata route unavailable; no file absence conclusion.'}|ConvertTo-Json -Depth 6
 exit 2
} finally {if($reader){$reader.Dispose()};if($command){$command.Dispose()};if($connection){$connection.Dispose()}}
