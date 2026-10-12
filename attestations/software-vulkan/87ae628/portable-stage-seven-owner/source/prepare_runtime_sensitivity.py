"""Select a source-derived owned node with no later owned replacement to a root."""
import hashlib, importlib.util, json, sys
from collections import deque
from pathlib import Path
sys.dont_write_bytecode = True
root = Path.cwd()
work = Path(__file__).resolve().parent
source = root/'src/sirius/kernels/retained_program.py'
spec = importlib.util.spec_from_file_location('producer',source)
module = importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
plans = module.build_stage_seven_owner_plan()
result = {}
for family,builder in (('general',module.build_endpoint_program),('schwarzschild',module.build_schwarzschild_endpoint_program)):
    p = builder(parallel=True,with_plan=True)
    plan = p['source_plan']; writes = plan['writes']; ops = plan['operations']
    slot,pair = max(enumerate(plans[family]['writes']),key=lambda item:item[1][1])
    ordinal = pair[1];node = writes[ordinal]
    live = set(writes); mapped = {writes[e] for _,e in plans[family]['writes']}
    def dependencies(op):
        code,a,b,c = op
        return () if code in (0,1) else (a,)+(() if code in (6,7,8) else (b,))+((c,) if code in (10,11) else ())
    parents = {n:[] for n in writes}
    for n in writes:
        for dep in dependencies(ops[n]):parents[dep].append(n)
    paths = {node:[node]};queue=deque([node])
    while queue:
        n=queue.popleft()
        for parent in parents[n]:
            assert parent not in mapped, 'highest mapped ordinal has a mapped consumer'
            if parent not in paths:paths[parent]=paths[n]+[parent];queue.append(parent)
    roots = [(i,n) for i,n in enumerate(plan['outputs']) if n in paths]
    assert roots
    result[family]={'slot':slot,'endpoint_ordinal':ordinal,'source_node':node,'source_opcode':ops[node][0],
                    'reachable_roots':{str(i):paths[n] for i,n in roots},'mapped_replacement_on_path':False}
header=['// Source-derived highest owned ordinal; downstream paths contain no replayed node.', '#pragma once']
for family,item in result.items():header.append(f'inline constexpr std::uint32_t k{family.title()}SensitivitySlot = {item["slot"]}u;')
(work/'runtime_sensitivity.h').write_text('\n'.join(header)+'\n')
(work/'runtime-sensitivity.json').write_text(json.dumps({'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'families':result,'scope':'Structural live-root reachability only. Actual conditional consumption and sensitivity require the controlled Vulkan execution.'},indent=2)+'\n')
print(json.dumps(result,indent=2))
