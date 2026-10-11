"""Restore pinned diagnostic source for the next finite saved-byte ISA join."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,subprocess
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def git(*a):return subprocess.check_output(['git',*a])
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
head=git('rev-parse','HEAD').decode().strip();assert head=='6be3132b7023f9d6497ab8f67dc779cf8eb4a80c'
assert git('rev-parse','@{u}').decode().strip()==head and not git('status','--porcelain').strip()
assert git('rev-parse','HEAD^{tree}')==git('rev-parse','b81e061^{tree}')
clock='0e2c7d64886b8218e9727057c067755b9089677f';views='30f2237b384ec20655db0539014a0fa2430f829c';restored={}
for revision,names in [(clock,('decode_saved_sqtt.py','run_cpu_decoder.py','rocprof_trace_decoder.h','trace_decoder_types.h','decoder-LICENSE','decoder-preparation.json','cpu-decoder-preparation.json','sdk-decode.cpp','aql-threadtrace.cpp','spm-index-result.json','clock-alignment-result.json')),(views,('trace-code-join.json','counter-findings.json'))]:
 for name in names:
  data=git('show',revision+':'+name);p=W/name;assert not p.exists();p.write_bytes(data)
  restored[name]={'revision':revision,'blob':git('rev-parse',revision+':'+name).decode().strip(),**seal(p)}
manifest=ROOT/'attestations/native-vulkan/6be3132/native-counter-clock/manifest.json';m=json.loads(manifest.read_text())
previous=next(v for v in m['records'] if v['origin']=='out/native-counter-clock/decoder-result.json');prior_decoder=ROOT/previous['recovery']['path'];assert seal(prior_decoder)=={k:previous[k] for k in ('bytes','sha256')}
trace=ROOT/'attestations/native-vulkan/6be3132/native-spm-cost/payloads/ee233c821d2e9692b800e17095ba71447e909f80f62180bcd66ac0dceccf05a9';assert seal(trace)=={'bytes':9581418,'sha256':trace.name}
library=ROOT/json.loads((W/'cpu-decoder-preparation.json').read_text())['extracted_library']['path'];assert seal(library)=={'bytes':199792,'sha256':'9aae97c42c748c095f2241421f9586a780d42b9d97678e2e061cf00d17cb150f'}
plan={'source_revision':head,'previous_goal_turn':'progress: actual offline decode/arithmetic evidence changed next action; useful source pushed, issue44 updated, scratch removed. No live handle carried forward.','trace':{'path':str(trace.relative_to(ROOT)),**seal(trace)},'prior_decoder_result':{'path':str(prior_decoder.relative_to(ROOT)),**seal(prior_decoder)},'prior_manifest':{'path':str(manifest.relative_to(ROOT)),**seal(manifest)},'restored_sources':restored,'limits':{'wall_seconds':45,'address_space_bytes':512*1024**2},'purpose':'Join exact saved code-load addresses/ELF/printed ISA to a valid callback; then inspect actual event-origin/clock records. No GPU acquisition or production change.'}
p=W/'plan.json';assert not p.exists();p.write_text(json.dumps(plan,indent=2)+'\n');print(json.dumps({'restored_sources':len(restored),'head':head,'library':seal(library),'previous_turn':'progress','native_launches':0}))
