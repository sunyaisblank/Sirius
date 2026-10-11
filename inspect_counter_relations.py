"""Compare complete observed raw sums and derived counter relationships."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,math,struct,subprocess
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def digest(values):return hashlib.sha256(struct.pack('<'+str(len(values))+'d',*values)).hexdigest()
plan=read(W/'plan.json');head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
assert head==plan['source_revision'] and not subprocess.check_output(['git','status','--porcelain'],text=True).strip()
trace=ROOT/plan['trace']['path'];assert seal(trace)=={k:plan['trace'][k] for k in ('bytes','sha256')}
paths=[Path(__file__),W/'plan.json',trace,W/'spm-index-result.json',W/'raw-spm-result.json',W/'derived-spm-result.json',W/'rdp-derived-header.h',W/'pal-gpu-perf-header.h']
inputs={str(p.relative_to(ROOT)):seal(p) for p in paths}
for name in ('spm-index-result.json','raw-spm-result.json','derived-spm-result.json','rdp-derived-header.h','pal-gpu-perf-header.h'):
 assert seal(W/name)=={k:plan['restored_sources'][name][k] for k in ('bytes','sha256')}
index=read(W/'spm-index-result.json');raw=read(W/'raw-spm-result.json');metadata=read(W/'derived-spm-result.json')
assert index['trace']==raw['trace']==metadata['trace']==seal(trace)
b=trace.read_bytes();n=raw['session']['timestamps'];assert n==11911
def data(c):
 assert c['compression']==0
 d=b[c['data_offset']:c['data_offset']+c['data_bytes']]
 assert len(d)==c['data_bytes'] and hashlib.sha256(d).hexdigest()==c['data_sha256']
 return d
sums={};instances={}
for c in index['chunks']:
 if c['name']!='SpmCounterData':continue
 assert c['version']==2 and c['header_bytes']==20
 pci,block,instance,event,width=struct.unpack_from('<5I',b,c['header_offset'])
 assert pci==raw['session']['pci_id'] and width==2
 values=struct.unpack('<'+str(n)+'H',data(c));key=(block,event)
 assert instance not in instances.setdefault(key,[]);instances[key].append(instance)
 acc=sums.setdefault(key,[0]*n)
 for i,v in enumerate(values):acc[i]+=v
assert len(sums)==27 and sum(len(v) for v in instances.values())==183
sum_hashes={key:digest(values) for key,values in sums.items()}
derived={};identities={}
for counter in metadata['counters']:
 c=index['chunks'][counter['chunk_ordinal']];assert c['name']=='DerivedSpmCtr' and c['version']==1
 d=data(c);idx=counter['counter_index'];values=struct.unpack_from('<'+str(n)+'d',d)
 assert all(math.isfinite(v) for v in values)
 assert hashlib.sha256(d[:n*8]).hexdigest()==counter['sample_data_sha256']
 derived[idx]=values;identities[idx]=counter
assert set(derived)==set(range(33))
joins=[]
for idx,values in derived.items():
 h=digest(values);matches=[{'gpu_block':k[0],'event_index':k[1],'instances':sorted(instances[k])} for k,v in sum_hashes.items() if v==h]
 joins.append({'derived_index':idx,'display_name':identities[idx]['display_name'],'usage_type_encoding':identities[idx]['usage_type_encoding'],'complete_raw_sum_byte_matches':matches,'scope':'Complete observed raw-instance sum converted to binary64 equals stored derived values. Multiple all-zero matches are ambiguous; equality alone is not a hardware counter name/formula/filter proof.'})
ratios=[]
for idx,num,den in ((0,2,1),(4,6,5),(8,10,9),(14,16,15),(18,20,19),(23,25,24),(26,28,24),(27,29,24)):
 assert num in identities[idx]['component_indices'] and den in identities[idx]['component_indices']
 assert identities[idx]['usage_type_encoding']==1
 expected=[];nonzero=0;errors=[];bit_differences=[]
 for i,(numerator,denominator,stored) in enumerate(zip(derived[num],derived[den],derived[idx])):
  value=(numerator/denominator)*100 if denominator!=0 else 0.0
  if denominator!=0:nonzero+=1
  expected.append(value);errors.append(abs(value-stored))
  if struct.pack('<d',value)!=struct.pack('<d',stored):bit_differences.append(i)
 ratios.append({'index':idx,'name':identities[idx]['display_name'],'numerator_index':num,'denominator_index':den,'observed_relation':'(numerator / denominator) * 100; zero denominator -> zero','nonzero_denominator_samples':nonzero,'maximum_absolute_error':max(errors),'binary64_different_samples':len(bit_differences),'first_different_samples':bit_differences[:8]})
result={'source_revision':head,'input_seals':inputs,'trace':seal(trace),'samples_per_stream':n,'raw_block_event_groups':len(sums),'raw_to_derived_observed_joins':joins,'ratio_observed_joins':ratios,'scope':'Finite whole-stream equality and arithmetic checks on this immutable trace only. No hardware/formula/filter universality, clock/event interval mapping, exclusive ownership, removable bottleneck, ordinary performance gain or scientific/release qualification.'}
for p,v in inputs.items():assert seal(ROOT/p)==v
p=W/'counter-relations.json';assert not p.exists();p.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({'complete_raw_sum_matches':[{'index':v['derived_index'],'name':v['display_name'],'matches':v['complete_raw_sum_byte_matches']} for v in joins if v['complete_raw_sum_byte_matches']],'ratios':ratios}))
