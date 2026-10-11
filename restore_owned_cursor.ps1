param([Parameter(Mandatory=$true)][string]$Output,[Parameter(Mandatory=$true)][string]$CursorSource,[Parameter(Mandatory=$true)][string]$CursorHash)
$ErrorActionPreference='Stop'
$receipt=[ordered]@{pass=$false;scope='Restore only this recorded cursor token, preserving any externally changed position; no clicks, input, activation or desktop settings.'}
try {
 $p=Join-Path $Output 'cursor-state.json'
 if (!(Test-Path -LiteralPath $p)) {$receipt.pass=$true;$receipt.result='no-cursor-token'}
 else {
  if ((Get-FileHash -LiteralPath $CursorSource -Algorithm SHA256).Hash.ToLowerInvariant() -ne $CursorHash) {throw 'Cursor helper source mismatch'}
  Add-Type -TypeDefinition ([IO.File]::ReadAllText($CursorSource))
  $s=Get-Content -Raw -LiteralPath $p | ConvertFrom-Json
  if ($s.before.Count -ne 2 -or $s.last.Count -ne 2 -or $s.target.Count -ne 2) {throw 'Cursor token dimensions invalid'}
  $receipt.result=[SiriusCursor]::Restore([int]$s.before[0],[int]$s.before[1],[int]$s.last[0],[int]$s.last[1],[int]$s.target[0],[int]$s.target[1])
  $receipt.pass=$receipt.result -in @('restored','already-original','external-position-preserved')
 }
} catch {$receipt.error=$_.Exception.Message}
try {$receipt | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $Output 'cursor-recovery.json') -Encoding UTF8}
catch {[Console]::Error.WriteLine('Cursor recovery receipt failed: '+$_.Exception.Message)}
# The caller always continues into owned process cleanup, then rejects a false
# receipt. A failed cursor restoration must not bypass tool termination.
