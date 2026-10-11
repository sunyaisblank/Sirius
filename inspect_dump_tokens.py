"""Read every retained debug callback byte; report grammar before interpreting it."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import collections,hashlib,json,struct
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
r=json.loads((W/'dump-result.json').read_text());inputs={str(p.relative_to(ROOT)):seal(p) for p in (Path(__file__),W/'dump-result.json',W/'dump-bootstrap.json',W/'xgl-marker-header.h')}
types={};counts=collections.Counter();count=0;total=0;backwards=0;last=None;all_digest=hashlib.sha256()
for entry in r['token_files']:
 p=ROOT/entry['path'];assert seal(p)=={k:entry[k] for k in ('bytes','sha256')};inputs[entry['path']]=seal(p);b=p.read_bytes();offset=0;all_digest.update(b);total+=len(b)
 while offset<len(b):
  assert offset+16<=len(b);clock,n,m=struct.unpack_from('<qII',b,offset);offset+=16;assert n<=4096 and m<=16384 and offset+n+m<=len(b)
  kind=b[offset:offset+n];offset+=n;info=b[offset:offset+m];offset+=m;count+=1;counts[kind.hex()]+=1
  if last is not None and clock<last:backwards+=1
  last=clock
  value=types.setdefault(kind.hex(),{'type_utf8':kind.decode('utf8','replace'),'first':[],'last':None,'minimum_time':clock,'maximum_time':clock})
  item={'ordinal':count-1,'time':clock,'info_utf8':info.decode('utf8','replace'),'info_hex':info.hex()}
  if len(value['first'])<6:value['first'].append(item)
  value['last']=item;value['minimum_time']=min(value['minimum_time'],clock);value['maximum_time']=max(value['maximum_time'],clock)
 assert offset==len(b)
assert count==r['callback_count'] and total==r['total_token_bytes'] and dict(counts)==r['type_counts_hex']
result={'source_revision':r['source_revision'],'input_seals':inputs,'records':count,'bytes':total,'concatenated_binary_sha256':all_digest.hexdigest(),'backward_time_records':backwards,'types':types,'scope':'Complete retained callback-byte inspection. No undocumented grammar accepted, PAL marker reconstruction, event or display-clock origin, cost or scientific/release verdict.'}
p=W/'dump-token-inspection.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'records':count,'bytes':total,'backwards':backwards,'types':{v['type_utf8']:counts[k] for k,v in types.items()}}))
