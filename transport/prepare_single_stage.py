"""Prepare one fixed diagnostic pipeline per fresh process under unchanged guards."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import ast, hashlib, json
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
def read(path): return json.loads(path.read_text(encoding='utf-8-sig'))
def seal(path): return {'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
def one(text, old, new):
    assert text.count(old) == 1, (old, text.count(old))
    return text.replace(old, new)
def save(path, content):
    assert not path.exists(), path
    path.write_text(content)
assert read(WORK/'placement-query-bridge.json')['accepted'] is False
assert read(WORK/'retirement-run.json')['accepted'] and read(WORK/'native-retirement.json')['pass']
base_plan = read(WORK/'plan.json')
base_inputs = [Path(__file__), WORK/'native_clock_placement.cs', WORK/'placement-owner.ps1',
               WORK/'placement-child.ps1', WORK/'run_placement_query.py', WORK/'source-review.json',
               WORK/'placement-query-bridge.json', WORK/'placement-query-emergency.json',
               WORK/'retirement-run.json', WORK/'native-retirement.json']
for index, stage in enumerate(('Transport', 'Endpoint')):
    target = WORK/stage.lower()
    target.mkdir()
    plan = dict(base_plan)
    plan['input_seals'] = dict(base_plan['input_seals'])
    plan['input_seals'].update({str(p.relative_to(ROOT)): seal(p) for p in base_inputs})
    plan['selected_stage'] = stage
    plan['purpose'] = ('One fixed '+stage+' diagnostic pipeline in a fresh native process. '
                       'Unchanged 50s sampled512MiB child guard/45s intrinsic query watchdog/120s bridge. '
                       'Combined query was invalidated during Endpoint compilation by RSS; fresh residency is a different query design. '
                       'No dispatch, timing benefit, numerical, scientific or release verdict.')
    for name in ('retained_transport.spv','retained_endpoint.spv','emergency_native_cleanup.ps1',
                 'run_native_retirement.py','verify_native_retirement.ps1'):
        p = target/name
        p.write_bytes((WORK/name).read_bytes())
        plan['input_seals'][str(p.relative_to(ROOT))] = seal(p)
    save(target/'plan.json', json.dumps(plan, indent=2)+'\n')
    source = (WORK/'native_clock_placement.cs').read_text()
    source = one(source, 'for (int i = 0; i < paths.Length; ++i)',
                 'for (int i = '+str(index)+'; i < '+str(index+1)+'; ++i)')
    source = one(source, '"pipelines_created", 2', '"pipelines_created", 1')
    source = one(source, 'Exact unadopted native-default phase-observer modules;',
                 'One fixed '+stage+' unadopted native-default phase-observer module in a fresh process;')
    save(target/'native_clock_placement.cs', source)
    owner = (WORK/'placement-owner.ps1').read_text()
    owner = one(owner, '$result.pipelines_created -ne 2', '$result.pipelines_created -ne 1')
    owner = one(owner, 'query for two exact diagnostic modules', 'query for one fixed '+stage+' diagnostic module in a fresh process')
    save(target/'placement-owner.ps1', owner)
    save(target/'placement-child.ps1', (WORK/'placement-child.ps1').read_text())
    bridge = (WORK/'run_placement_query.py').read_text()
    bridge = one(bridge, "WORK/'prepare_placement.py'", "WORK/'prepare_single_stage.py'")
    bridge = one(bridge, "result['pipelines_created']==2", "result['pipelines_created']==1")
    bridge = one(bridge, "len(result['observations'])==2", "len(result['observations'])==1 and result['observations'][0]['stage']==plan['selected_stage']")
    bridge = one(bridge, "'pipeline_count':2", "'pipeline_count':1")
    bridge = one(bridge, "assert result['enabled_extensions']==", "assert result['extension_count']==len(result['extensions'])\n  assert result['enabled_extensions']==")
    save(target/'run_placement_query.py', bridge)
    ast.parse(bridge)
    (target/'prepare_single_stage.py').write_bytes(Path(__file__).read_bytes())
    outputs = {name:seal(target/name) for name in ('plan.json','native_clock_placement.cs',
               'placement-owner.ps1','placement-child.ps1','run_placement_query.py',
               'emergency_native_cleanup.ps1','prepare_single_stage.py')}
    save(target/'preparation.json', json.dumps({'selected_stage':stage,
         'base_input_seals':{str(p.relative_to(ROOT)):seal(p) for p in base_inputs},
         'outputs':outputs,'scope':'Source preparation only; no helper/compiler/native invocation.'},indent=2)+'\n')
    print(json.dumps({'stage':stage,'outputs':outputs}))
