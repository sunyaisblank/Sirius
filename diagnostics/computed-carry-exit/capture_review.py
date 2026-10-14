"""Root-only explicit later action: unchanged canary collection then byte review."""
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
work = Path(__file__).resolve().parent
folder = Path(sys.argv[1]).resolve()
assert folder.parent == work and folder.name == 'capture-candidate-transport-normal'
subprocess.run([sys.executable, '-B', str(work / 'run_capture.py'),
                str(folder), sys.argv[2]], check=True)
subprocess.run([sys.executable, '-B', str(work / 'analyze_capture.py'), folder.name], check=True)
