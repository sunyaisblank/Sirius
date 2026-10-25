"""Observe the failed C3 prefix once to choose arithmetic/JIT/ownership repair."""
from pathlib import Path
import json
import subprocess

R=Path.cwd();W=Path(__file__).resolve().parent
assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert not (W/'actions.json').exists() and not (W/'capture').exists()
source=R/'out/native-cpu-scalar-zero-cohort/controls/runtime'
failed=json.loads((source/'control-3/owner.json').read_text())
assert failed['status']=='terminal' and failed['stop_reason']=='diagnostic_timeout'
assert not failed['remaining_owned_processes'] and not failed['cleanup_errors']
action=json.loads((source/'actions.json').read_text())['control-3']
inputs=set(action['frozen_inputs'])|{str(path.relative_to(R)) for path in W.iterdir() if path.is_file()}
inputs|={'/usr/bin/gdb',str((source/'control-3/owner.json').relative_to(R)),
    'attestations/qualification/fb50bb2a/native-cpu-scalar-zero-control3-timeout/MANIFEST.json'}
action.update(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),kind='diagnostic',
    limits={'seconds':60,'rss_kib':4194304,'reserve_kib':2097152},expected_exit=0,
    argv=['/usr/bin/python3','-B',str(W/'capture.py'),'capture','{folder}'],frozen_inputs=sorted(inputs),
    scope='One diagnostic original failed-C3 prefix with a fixed15s stack/first scalar input capture, then deliberate termination. Decision: distinguish native arithmetic from LLVM/JIT or ownership stall before changing implementation. Stricter60s outer diagnostic bound, same4GiB/2GiB reserve; no original control retry/acceptance, performance comparison or frame.')
(W/'actions.json').write_text(json.dumps({'capture':action},indent=2)+'\n')
print(json.dumps({'diagnostic_prepared':True,'fixed_capture_seconds':15,'executed':False}))
