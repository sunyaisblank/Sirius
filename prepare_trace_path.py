"""Give the protected payload its original .rgp file type for viewer opening."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,shutil
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
plan_path=WORK/'counter-view-plan.json';bridge_path=WORK/'run_profiler_counter_inventory.py'
plan=json.loads(plan_path.read_text());protected=dict(plan['trace']);origin=ROOT/protected['path']
assert seal(origin)=={k:protected[k] for k in ('bytes','sha256')}
trace=WORK/'counter-enabled-interval.rgp';assert not trace.exists();shutil.copy2(origin,trace);assert seal(trace)==seal(origin)
plan['protected_trace']=protected;plan['trace']={'path':str(trace.relative_to(ROOT)),**seal(trace)}
initial={'plan':{'bytes':plan_path.stat().st_size,'sha256':seal(plan_path)['sha256']},'bridge':seal(bridge_path),'scope':'Before adapting extensionless protected-payload reference to a task-local .rgp viewer copy; no viewer had been launched.'}
plan['pre_viewer_file_type_adaptation']=initial
bridge=bridge_path.read_text();old="paths += [trace,WORK/'counter-view-plan.json',WORK/'prepare_counter_views.py',ROOT/plan['prior_manifest']['path']]"
new="paths += [trace,ROOT/plan['protected_trace']['path'],WORK/'counter-view-plan.json',WORK/'prepare_counter_views.py',WORK/'prepare_trace_path.py',ROOT/plan['prior_manifest']['path']]"
assert bridge.count(old)==1
bridge_path.write_text(bridge.replace(old,new));plan_path.write_text(json.dumps(plan,indent=2)+'\n')
print(json.dumps({'viewer_trace':str(trace.relative_to(ROOT)),'protected_payload_unchanged':True,'trace':seal(trace)}))
