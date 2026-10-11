"""Use the pre-existing ordinary static-query budget for Endpoint compilation only."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import ast, hashlib, json
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
BASE=WORK/'endpoint';TARGET=WORK/'endpoint-query-budget'
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p): return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def one(text,old,new):
    assert text.count(old)==1,(old,text.count(old))
    return text.replace(old,new)
def save(name,text): p=TARGET/name;assert not p.exists();p.write_text(text)
assert not read(BASE/'placement-query-bridge.json')['accepted']
assert read(BASE/'retirement-run.json')['accepted'] and read(BASE/'native-retirement.json')['pass']
previous=ROOT/'attestations/native-vulkan/3045067/ordinary-pipeline-generated-code/run-query/owner.json'
endpoint=ROOT/'attestations/native-vulkan/3045067/ordinary-pipeline-generated-code/run-endpoint-query/owner.json'
for p in (previous,endpoint):
    r=read(p);assert r['status']=='completed' and r['outer_guard_seconds']==90 and r['rss_guard_bytes']==4294967296
assert read(previous)['peak_sampled_rss_bytes']>536870912
TARGET.mkdir()
plan=read(BASE/'plan.json')
plan['input_seals']=dict(plan['input_seals'])
inputs=[Path(__file__),previous,endpoint,BASE/'native_clock_placement.cs',BASE/'placement-owner.ps1',
        BASE/'placement-child.ps1',BASE/'run_placement_query.py',BASE/'source-review.json',
        BASE/'placement-query-bridge.json',BASE/'placement-query-emergency.json',
        BASE/'retirement-run.json',BASE/'native-retirement.json']
plan['input_seals'].update({str(p.relative_to(ROOT)):seal(p) for p in inputs})
plan['limits']=dict(plan['limits'],native_seconds=90,native_rss_bytes=4294967296,native_intrinsic_watchdog_seconds=85)
plan['purpose']=('Endpoint-only ordinary static generated-code query under the pre-existing90s/4GiB diagnostic compiler budget; '
                 'added85s parent-birth-held intrinsic watchdog and unchanged120s bridge. Both512MiB stops remain invalidated. '
                 'Historical original module compilation sampled577490944B, already beyond512MiB. '
                 'This changes the tooling-query budget only; no GPU dispatch, scientific/runtime/performance acceptance or original rendering governor change.')
for name in ('retained_transport.spv','retained_endpoint.spv','emergency_native_cleanup.ps1',
             'run_native_retirement.py','verify_native_retirement.ps1','placement-child.ps1'):
    (TARGET/name).write_bytes((BASE/name).read_bytes())
    plan['input_seals'][str((TARGET/name).relative_to(ROOT))]=seal(TARGET/name)
save('plan.json',json.dumps(plan,indent=2)+'\n')
source=one((BASE/'native_clock_placement.cs').read_text(),'timer.Elapsed.TotalSeconds>45','timer.Elapsed.TotalSeconds>85')
save('native_clock_placement.cs',source)
owner=one((BASE/'placement-owner.ps1').read_text(),'[int]$Seconds=50,[long]$RssBytes=536870912','[int]$Seconds=90,[long]$RssBytes=4294967296')
owner=one(owner,'50s native polling/sample512MiB childRSS guard;120s bridge plus separate cleanup bounds.','Existing90s native polling/sample4GiB childRSS tooling-query guard;85s intrinsic query watchdog;120s bridge plus separate cleanup bounds. No rendering or numerical governor changed.')
save('placement-owner.ps1',owner)
bridge=(BASE/'run_placement_query.py').read_text()
bridge=one(bridge,"WORK/'prepare_single_stage.py'","WORK/'prepare_endpoint_query_budget.py'")
bridge=one(bridge,"'-Seconds','50','-RssBytes','536870912'","'-Seconds','90','-RssBytes','4294967296'")
bridge=one(bridge,"owner['outer_guard_seconds']==50 and owner['rss_guard_bytes']==536870912 and owner['peak_sampled_rss_bytes']<=536870912","owner['outer_guard_seconds']==90 and owner['rss_guard_bytes']==4294967296 and owner['peak_sampled_rss_bytes']<=4294967296")
save('run_placement_query.py',bridge);ast.parse(bridge)
(TARGET/'prepare_endpoint_query_budget.py').write_bytes(Path(__file__).read_bytes())
outputs={name:seal(TARGET/name) for name in ('plan.json','native_clock_placement.cs','placement-owner.ps1',
         'placement-child.ps1','run_placement_query.py','emergency_native_cleanup.ps1','prepare_endpoint_query_budget.py')}
save('preparation.json',json.dumps({'base_input_seals':{str(p.relative_to(ROOT)):seal(p) for p in inputs},
     'outputs':outputs,'scope':'Source preparation only; no helper/compiler/native invocation. Existing tooling-query budget; no acceptance-governor change.'},indent=2)+'\n')
print(json.dumps(outputs))
