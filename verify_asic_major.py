"""Check the recorded major version needed by the public normalization branch."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,struct
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
plan=json.loads((W/'plan.json').read_text());trace=ROOT/plan['trace']['path'];assert seal(trace)=={k:plan['trace'][k] for k in ('bytes','sha256')}
idx=json.loads((W/'spm-index-result.json').read_text());c=next(v for v in idx['chunks'] if v['name']=='AsicInfo');b=trace.read_bytes();d=b[c['data_offset']:c['data_offset']+c['data_bytes']];assert hashlib.sha256(d).hexdigest()==c['data_sha256']
assert struct.unpack_from('<3i',d,64)==(1,12,2)
major=struct.unpack_from('<H',d,104)[0];assert major==11
paths=[Path(__file__),trace,W/'plan.json',W/'spm-index-result.json',W/'rdp-file-format.h',W/'rdp-derived-database.cpp',W/'normalization-result.json']
result={'source_revision':plan['source_revision'],'input_seals':{str(p.relative_to(ROOT)):seal(p) for p in paths},'recorded_gfx_major':major,'byte_offset':104,'public_rdp_field_encoding':'TraceGfxIpLevel begins with uint16_t major. Only this major field is needed for isGfx10OrHigher.','gfx10_or_higher_branch_selected':True,'public_normalization_num_simds':48,'scope':'Major-version branch for source-defined counter normalization only. No assertion about other GFX version fields, physical CU inventory, packaged DLL equivalence, epoch or qualification.'}
(W/'normalization-gfx-major.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'gfx_major':major,'num_simds':48}))
