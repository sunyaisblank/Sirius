"""Run the four remaining original controls in order, once without diagnostics."""
from pathlib import Path
import subprocess

W = Path(__file__).resolve().parent
for i in range(4, 8):
    subprocess.run(['/usr/bin/python3', '-B', str(W / 'runtime-unobserved/run_owned.py'), 'control-' + str(i)], check=True)
subprocess.run(['/usr/bin/python3', '-B', str(W / 'readers-unobserved/read_controls.py')], check=True)
