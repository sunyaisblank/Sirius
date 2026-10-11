"""Reuse the owned profiler inventory on the exact preserved counter-enabled trace."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import ast,hashlib,json,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
OLD='evidence/issue-44-native-executed-cost-6be3132'
SPM='evidence/issue-44-native-spm-6be3132'
assert ROOT==Path('/home/astra/.project/Sirius') and WORK==ROOT/'out/native-counter-views'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def git(*a):return subprocess.check_output(['git',*a],text=True).strip()
def blob(tag,n):return subprocess.check_output(['git','show',tag+':'+n])
def replace(s,o,n):
 assert s.count(o)==1,(o[:100],s.count(o))
 return s.replace(o,n)
head=git('rev-parse','HEAD');assert head==git('rev-parse','@{u}') and not git('status','--porcelain')
assert git('rev-parse','HEAD^{tree}')==git('rev-parse','b81e061^{tree}')
prior=ROOT/'attestations/native-vulkan/6be3132/native-spm-cost/manifest.json'
m=read(prior);assert m['source_revision']==head and m['source_snapshot']['commit']==git('rev-parse',SPM)
v=next(v for v in m['records'] if v['origin'].endswith('/spm-original-legacy-interval.rgp'))
assert v['recovery']['kind']=='local-archive-payload'
trace=ROOT/v['recovery']['path'];assert seal(trace)=={k:v[k] for k in ('bytes','sha256')}
assert v['bytes']==9581418 and v['sha256']=='ee233c821d2e9692b800e17095ba71447e909f80f62180bcd66ac0dceccf05a9'
outputs={};ancestry={}
def save(n,b):
 assert n not in outputs and not (WORK/n).exists();outputs[n]=b.encode() if isinstance(b,str) else b
for n in ('tools-extraction.json','tools-gui-extraction.json','emergency_native_cleanup.ps1'):
 b=blob(SPM,n);save(n,b);ancestry[n]={'revision':git('rev-parse',SPM),'path':n,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
for n in ('profiler_interval_inventory.ps1','run_profiler_interval_inventory.py'):
 b=blob(OLD,n);ancestry[n]={'revision':git('rev-parse',OLD),'path':n,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
ps=blob(OLD,'profiler_interval_inventory.ps1').decode()
ps=ps.replace('profiler-interval-inventory','profiler-counter-inventory')
ps=ps.replace('Event timing','Wavefront occupancy')
save('profiler_counter_inventory.ps1',ps)
b=blob(OLD,'run_profiler_interval_inventory.py').decode()
b=b.replace('profiler-interval-inventory','profiler-counter-inventory')
b=b.replace('profiler_interval_inventory.ps1','profiler_counter_inventory.ps1')
b=b.replace('interval-inventory-source-review.json','counter-inventory-source-review.json')
b=b.replace('Event timing','Wavefront occupancy')
begin=b.index("trace=WORK/'native-capture-interval-console")
end=b.index("paths += [ROOT/v['path']",begin)
b=b[:begin]+'''plan=read(WORK/'counter-view-plan.json');assert plan['source_revision']==head
trace=ROOT/plan['trace']['path'];assert trace.stat().st_size==plan['trace']['bytes'] and sha(trace)==plan['trace']['sha256']
paths += [trace,WORK/'counter-view-plan.json',WORK/'prepare_counter_views.py',ROOT/plan['prior_manifest']['path']]
'''+b[end:]
save('run_profiler_counter_inventory.py',b)
source={'source_revision':head,'branch':git('branch','--show-current'),'issue':44,
        'previous_goal_turn_disposition':'progress: actual counters acquired and validated; complete evidence/source pushed and finished scratch removed.',
        'scope':'Saved counter-enabled trace cost-analysis prerequisite; no production/scientific source or stream change.'}
save('source.json',json.dumps(source,indent=2)+'\n')
plan={'source_revision':head,'prior_manifest':{'path':str(prior.relative_to(ROOT)),**seal(prior)},
      'trace':{'path':str(trace.relative_to(ROOT)),**seal(trace)},'ancestral_sources':ancestry,
      'nav_history':['EVENTS','Wavefront occupancy'],'original_science_unchanged':True,
      'guard_seconds':60,'bridge_seconds':120,'tool_rss_bytes':1073741824,
      'scope':'One owned signed RGP saved-trace viewer; exact prior data/control pedigree. Inventory actual Wavefront occupancy controls and capture bounded owned-window PNG. No selected-dispatch identity, counter interpretation, source/event/clock join or cost attribution is assumed. No RDS/RDP/scientific/capture execution.'}
save('counter-view-plan.json',json.dumps(plan,indent=2)+'\n')
for n,b in outputs.items():
 if n.endswith('.py'):ast.parse(b,filename=n)
for n,b in outputs.items():(WORK/n).write_bytes(b)
print(json.dumps({'prepared_files':len(outputs),'trace':plan['trace'],'navigation':plan['nav_history']}))
