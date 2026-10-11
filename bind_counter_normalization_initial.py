"""Bind the observed LDS divisor to recorded ASIC fields and pinned public consumers."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,re,struct,subprocess,urllib.request
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def gh(path):return json.loads(subprocess.check_output(['gh','api',path]))
rdp='d63073a854c8ccc4b10edbf27559bb698c5685b3';gpa='5b7a7544eac3ebc7f0f62ff878c69f8f893f3583'
origins=[('GPUOpen-Tools/radeon_developer_panel',rdp,'source/trace/inc/rgp_file_format.h','rdp-file-format.h','23f707aa03fbaf6e0b0fef9d2fe3bde6a178e71113366da0f743577847967a5e'),('GPUOpen-Tools/gpu_performance_api',gpa,'source/gpu_perf_api_counters/gpa_counter_context_handler.cc','gpa-context-handler.cc','107a44735d84c9700b9ddde99092668f9b6581dd7a2ef48e34e0071d269a2e9e'),('GPUOpen-Tools/gpu_performance_api',gpa,'source/gpu_perf_api_counter_generator/gpa_derived_counter_evaluator.hpp','gpa-counter-evaluator.hpp','de76ca75f45d7c019448ce3fe79bad6928ba6668346b969fcf9f7fecc9d6f8fa'),('GPUOpen-Tools/gpu_performance_api',gpa,'source/gpu_perf_api_common/gpa_hw_info.cc','gpa-hardware-info.cc','1628ffd499541d29cb904bef15ea810fb769da5e8915a69143d23a4d741ffc40'),('GPUOpen-Tools/gpu_performance_api',gpa,'LICENSE.txt','gpa-LICENSE',None),('GPUOpen-Tools/radeon_developer_panel',rdp,'LICENSE.txt','rdp-LICENSE',None)]
sources=[]
for repo,revision,origin,name,expected_sha in origins:
 info=gh('repos/'+repo+'/contents/'+origin+'?ref='+revision);assert info['type']=='file' and info['size']<65536
 url='https://raw.githubusercontent.com/'+repo+'/'+revision+'/'+origin
 with urllib.request.urlopen(url,timeout=30) as r:data=r.read(65537)
 assert len(data)==info['size'] and hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==info['sha']
 if expected_sha:assert hashlib.sha256(data).hexdigest()==expected_sha
 p=W/name;assert not p.exists();p.write_bytes(data)
 sources.append({'path':str(p.relative_to(ROOT)),'url':url,'repository':repo,'revision':revision,'repository_path':origin,'git_blob':info['sha'],**seal(p)})
plan=json.loads((W/'plan.json').read_text());trace=ROOT/plan['trace']['path'];assert seal(trace)=={k:plan['trace'][k] for k in ('bytes','sha256')}
index=json.loads((W/'spm-index-result.json').read_text());c=next(v for v in index['chunks'] if v['name']=='AsicInfo');assert c['version']==3 and c['data_bytes']==608 and c['data_sha256']=='a80c73cbe11275d7da7f7aed11ec23a3a9b74d39ea0284f4e87d268d02d1ec9c'
data=trace.read_bytes()[c['data_offset']:c['data_offset']+c['data_bytes']];assert hashlib.sha256(data).hexdigest()==c['data_sha256']
se,cu,simd=struct.unpack_from('<3i',data,64);assert (se,cu,simd)==(1,12,2)
public=(W/'rdp-derived-database.cpp').read_text();assert 'num_cu *= 2;' in public and 'num_simd                 = num_cu * asic_parameters.simd_per_compute_unit;' in public
assert 'kGpaHardwareAttributeNumSimds' in (W/'gpa-context-handler.cc').read_text() and 'SetNumberSimds' in (W/'gpa-context-handler.cc').read_text()
assert 'NUM_SIMDS' in (W/'gpa-counter-evaluator.hpp').read_text() and 'hw_info.GetNumberSimds()' in (W/'gpa-counter-evaluator.hpp').read_text()
num_simds=se*cu*2*simd;assert num_simds==48
paths=[Path(__file__),W/'plan.json',trace,W/'spm-index-result.json',W/'rdp-derived-database.cpp',W/'public-formula-corrected-result.json']+[ROOT/v['path'] for v in sources]
result={'source_revision':plan['source_revision'],'input_seals':{str(p.relative_to(ROOT)):seal(p) for p in paths},'primary_sources':sources,'asic_chunk':c,'recorded_asic_fields':{'shader_engines':se,'compute_units_per_shader_engine':cu,'simds_per_compute_unit':simd,'field_byte_offsets':[64,68,72]},'source_rule':'RDP doubles num_cu for gfx>=10, then multiplies simd_per_compute_unit; hands kGpaHardwareAttributeNumSimds to GPA. GPA context stores it and the derived evaluator reads NUM_SIMDS from that attribute.','num_simds':num_simds,'formula_receipt':seal(W/'public-formula-corrected-result.json'),'scope':'Recorded ASIC metadata joined to pinned public RDP/GPA normalization rule. Supersedes the earlier unbound-divisor limitation for that source rule only; does not attest packaged DLL build/formula/filter, physical inventory, event attribution, ordinary benefit or qualification.'}
(W/'normalization-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'num_simds':num_simds,'primary_sources':len(sources),'asic_fields':[se,cu,simd]}))
