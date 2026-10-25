"""Stop the original frame entry after first positive route coverage; no timing verdict."""
from pathlib import Path
import json
import subprocess

R=Path.cwd()
W=Path(__file__).resolve().parent
prepared=json.loads((W/'prepared.json').read_text())
code=subprocess.run(prepared['gdb_argv'],cwd=R).returncode
observation=json.loads((W/'observe/native-observation.json').read_text())
result={'debugger_exit':code,'capture_complete':observation['capture_complete'],
        'failures':observation['failures'],'frame_qualification':False,
        'scope':'One bounded first-route observation of the original frame test with unchanged scene, binaries, cache policy and memory/reserve limits. Inferior deliberately stops after first native3/SIMD3/dense4 coverage; no completed frame, timing comparison, numerical or release acceptance.'}
(W/'result.json').write_text(json.dumps(result,indent=2)+'\n')
raise SystemExit(0 if code==0 and result['capture_complete'] and not result['failures'] else 1)
