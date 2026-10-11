"""Restore exact previous diagnostic sources and protected results for one code-end check."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,subprocess
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def git(*args):return subprocess.check_output(['git',*args])
head=git('rev-parse','HEAD').decode().strip();assert head=='6be3132b7023f9d6497ab8f67dc779cf8eb4a80c' and git('rev-parse','@{u}').decode().strip()==head
assert not git('status','--porcelain').strip() and not any(p for p in W.iterdir() if p.name!='prepare_code_end.py')
assert git('rev-parse','HEAD^{tree}')==git('rev-parse','b81e061^{tree}')
revision='06e85e0630cef0bfa3ec9fcef4cb44c59a63db75';restored={}
names=('decode_with_normalized_branches.py','run_normalized_isa_decoder.py','rocprof_trace_decoder.h','trace_decoder_types.h','sdk-decode.cpp','sdk-disassembly.hpp','aql-threadtrace.cpp','cpu-decoder-preparation.json','decoder-preparation.json','decoder-LICENSE','spm-index-result.json','trace-code-join.json','code-end-byte-result.json','rocm-gfx11-sopp-test.txt','rocm-llvm-LICENSE.txt','pal-code-load-header.h','pal-code-load-producer.cpp','pal-code-uploader.cpp')
for name in names:
 p=W/name;data=git('show',revision+':'+name);assert not p.exists();p.write_bytes(data);restored[name]={'revision':revision,'blob':git('rev-parse',revision+':'+name).decode().strip(),**seal(p)}
D=ROOT/'attestations/native-vulkan/6be3132/native-sqtt-source-clock';manifest=D/'manifest.json';assert seal(manifest)=={'bytes':155076,'sha256':'08e33d89fb29db7c84a3f5f0b43c4b568d9db49fe05433fa365ad1daa6d08a75'}
m=json.loads(manifest.read_text());prior={}
for name in ('saved-isa2-map.json','isa2-result.json','isa2-run.json'):
 v=next(v for v in m['records'] if v['origin']=='out/native-sqtt-source-clock/'+name);p=ROOT/v['recovery']['path'];assert seal(p)=={k:v[k] for k in ('bytes','sha256')};prior[name]={'path':str(p.relative_to(ROOT)),**seal(p)}
trace=ROOT/'attestations/native-vulkan/6be3132/native-spm-cost/payloads/ee233c821d2e9692b800e17095ba71447e909f80f62180bcd66ac0dceccf05a9';assert seal(trace)=={'bytes':9581418,'sha256':trace.name}
library=ROOT/json.loads((W/'cpu-decoder-preparation.json').read_text())['extracted_library']['path'];assert seal(library)=={'bytes':199792,'sha256':'9aae97c42c748c095f2241421f9586a780d42b9d97678e2e061cf00d17cb150f'}
plan={'source_revision':head,'previous_turn':'progress: actual15977-PC first-wave binding, reviewed preserved/pushed sources, issue44 updated and completed scratch removed. No live handle carried forward.','trace':{'path':str(trace.relative_to(ROOT)),**seal(trace)},'prior_manifest':{'path':str(manifest.relative_to(ROOT)),**seal(manifest)},'prior':prior,'restored_sources':restored,'limits':{'wall_seconds':45,'address_space_bytes':512*1024**2},'purpose':'Bind only the actually requested Transport entry-end padding from whole-byte source-bound allocated text; inspect actual remaining stitching and existing artifact recovery without GPU work.'}
(W/'plan.json').write_text(json.dumps(plan,indent=2)+'\n');print(json.dumps({'restored_files':len(restored),'head':head,'prior_map':prior['saved-isa2-map.json'],'previous_turn':'progress','native_launches':0}))
