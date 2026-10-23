"""Reuse the authenticated existing incremental builder."""
from pathlib import Path
import hashlib
_basis=Path.cwd()/'out/lavapipe-owner-bulk-zero/prepare_build.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest()=='4dc31d0875ad372df65a8604c69b3f8daaa08d57959ad3787d2cfcd69694c674'
exec(compile(_basis.read_text(),str(_basis),'exec'))
