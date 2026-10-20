from pathlib import Path
import json,subprocess
W=Path(__file__).resolve().parent
for row in json.loads((W/'build-commands.json').read_text()):
    subprocess.run(row['argv'],cwd=row['cwd'],check=True)
