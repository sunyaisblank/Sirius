"""Match public decoder reference values to chronological saved-payload bit patterns."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
plan=json.loads((W/'plan.json').read_text());trace=ROOT/plan['trace']['path'];assert seal(trace)=={k:plan['trace'][k] for k in ('bytes','sha256')}
idx=json.loads((W/'spm-index-result.json').read_text());c=next(v for v in idx['chunks'] if v['name']=='SqttData');b=trace.read_bytes();raw=b[c['data_offset']:c['data_offset']+c['data_bytes']]
assert hashlib.sha256(raw).hexdigest()==c['data_sha256'];d=json.loads((W/'decoder-result.json').read_text());pairs=d['realtime_records'];assert len(pairs)==11912
nibbles=''.join(format(v&15,'x')+format(v>>4,'x') for v in raw);start=0;positions=[];following=[]
for v in pairs:
 rt=v['realtime_clock'];assert 0<=rt<1<<36
 pattern='102'+format(rt,'09x')[::-1]
 pos=nibbles.find(pattern,start);assert pos>=0 and pos+12<len(nibbles)
 positions.append(pos*4+12);following.append(int(nibbles[pos+12],16));start=pos+12
assert positions[:2]==[3084,6324] and positions[-2:]==[6246996,6247252]
paths=[Path(__file__),W/'plan.json',trace,W/'spm-index-result.json',W/'decoder-result.json'];inputs={str(p.relative_to(ROOT)):seal(p) for p in paths}
result={'source_revision':plan['source_revision'],'input_seals':inputs,'payload':{'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()},'records_matched':len(positions),'observed_pattern':'12 LSB-first bits with value 0x201, then 36 LSB-first bits equal to returned realtime value; chronological nibble-aligned exact payload matches.','first_pattern_bit_positions':positions[:4],'last_pattern_bit_positions':positions[-4:],'next_nibble_nonzero_records':sum(v!=0 for v in following),'pattern_positions_sha256':hashlib.sha256(json.dumps(positions,separators=(',',':')).encode()).hexdigest(),'scope':'Finite raw-payload bit correspondence only; no decoder-binary analysis, published hardware token-layout proof, exact SPM sample identity, certified epoch/frequency, RGP event origin, exclusive cost or qualification.'}
(W/'clock-bit-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ('records_matched','first_pattern_bit_positions','last_pattern_bit_positions','next_nibble_nonzero_records')}))
