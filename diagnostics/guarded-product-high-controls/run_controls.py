"""Execute each unchanged original filter once, in dependency order."""
from pathlib import Path
import json, subprocess, sys

W = Path(__file__).resolve().parent
for index in range(1, 9):
    label = 'control-' + str(index)
    result = subprocess.run([sys.executable, '-B', str(W / 'run_owned.py'), label])
    assert result.returncode == 0, (label, result.returncode)
    owner = json.loads((W / label / 'owner.json').read_text())
    assert owner['status'] == 'terminal' and owner['passed']
