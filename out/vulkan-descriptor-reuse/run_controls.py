from pathlib import Path
import json, subprocess, sys

W = Path(__file__).resolve().parent
for i in range(2, 12):
    label = 'control-' + str(i)
    result = subprocess.run([sys.executable, '-B', str(W / 'run_owned.py'), label])
    assert result.returncode == 0, (label, result.returncode)
    owner = json.loads((W / label / 'owner.json').read_text())
    assert owner['passed'] and owner['status'] == 'terminal'
