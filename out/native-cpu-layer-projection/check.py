"""Verify actual local promotion; reject the measured interface cost."""
from pathlib import Path
import json,re

R=Path.cwd();W=Path(__file__).resolve().parent
d=json.loads((W/'cleanup-result.json').read_text());n=json.loads((W/'native-result.json').read_text())
assert all(x['local_allocas']==0 and x['masked_calls']==0 for x in d['projected_cores'])
text=(R/'bin/linux-gcc/native-cpu-layer-projection/dp.cleanup.ll').read_text()
cores=list(re.finditer(r'^define [^\n]*@(sirius_projected_layer_\d+)\([^\n]*\{\n(.*?)^}',text,re.M|re.S))
assert len(cores)==41
assert not any('ptr poison' in m[2] or re.search(r'ptr[^,]* %kernelContext_\d+(?:[, )])',m[2]) for m in cores)
record=dict(status='LOCAL_SSA_GATE_PASS_INTERFACE_COST_GATE_FAIL_STOP',
 local_projection=dict(cores=41,local_allocas_after_cleanup=0,masked_memory_calls_after_cleanup=0,poison_or_old_context_pointers=False,
 ordinary_core_loads=sum(x['ordinary_loads'] for x in d['projected_cores']),ordinary_core_stores=sum(x['ordinary_stores'] for x in d['projected_cores'])),
 strict_math_calls=sum(d['core_calls'].values()),strict_FP_inventory_unchanged=d['constrained_calls']==json.loads((W/'lower-result.json').read_text())['constrained_calls'],
 native_instructions_before=225609,native_instructions_after=n['instructions'],instruction_difference=n['instructions']-225609,
 continue_T_or_E_or_runtime=False,
 reason='Per-slot locals becomeSSA, but snapshot/copyback through41 compact buffers costs more native instructions. Do not run a new frame. Remove generic metadata/root escapes in a new admitted-only canonical DP copy, inline raw helpers without embedding shared math bodies, and normalize the single fully initialized owned context instead of copying it between layer boundaries.',
 scope='One bounded compiler prototype; no provider link/loading/numerical/speed/scientific/hardware/release qualification. Source semantic review covers initialized inactive slots, exact masked merges, readonly external views and all return copybacks under the pinned positive private-job contract.')
assert record['strict_math_calls']==1581 and record['strict_FP_inventory_unchanged']
p=W/'structural-gate.json';assert not p.exists();p.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))
