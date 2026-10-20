from pathlib import Path
import json,subprocess
W=Path(__file__).resolve().parent
for row in json.loads((W/'baseline-commands.json').read_text()):subprocess.run(row['argv'],cwd=row['cwd'],check=True)
subprocess.run([str(Path.cwd()/'bin/linux-gcc/vulkan-int64-admission-fix/baseline')],check=True)
