"""Join complete saved query text to ELF instruction bytes; never executes shaders."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,re,struct
WORK=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def elf_text(data):
    assert data[:6]==b'\x7fELF\x02\x01'
    offset=struct.unpack_from('<Q',data,40)[0]
    stride,count,names=struct.unpack_from('<HHH',data,58)
    assert stride==64 and names<count
    headers=[struct.unpack_from('<IIQQQQIIQQ',data,offset+i*stride) for i in range(count)]
    h=headers[names];strings=data[h[4]:h[4]+h[5]]
    found=[]
    for h in headers:
        end=strings.index(0,h[0]);name=strings[h[0]:end].decode('ascii')
        if name=='.text':found.append(data[h[4]:h[4]+h[5]])
    assert len(found)==1
    return found[0]
result={'source_revision':read(WORK/'placement/plan.json')['source_revision'],'stages':{},
        'scope':'Static ordinary generated-code placement only; no dispatch, clock execution/rate/conversion, numerical/scientific/release qualification or ordinary timing benefit. Sidecar observer perturbs intervals.'}
for folder,stage,site_count in [('placement','Transport',7)]:
    w=WORK/folder
    receipt=read(w/'clock-placement.json');bridge=read(w/'placement-query-bridge.json');owner=read(w/'run/owner.json')
    assert receipt['completed'] and bridge['accepted'] and owner['status']=='completed'
    assert read(w/'retirement-run.json')['accepted'] and read(w/'native-retirement.json')['pass']
    assert receipt['pipelines_created']==1 and receipt['dispatches']==0 and receipt['descriptor_bindings']==[0,1,2]
    observation=receipt['observations'][0];assert observation['stage']==stage
    module=w/('retained_'+stage.lower()+'.spv')
    assert seal(module)=={'bytes':observation['module_bytes'],'sha256':observation['module_sha256']}
    payloads={};uncaptured=[]
    for q in observation['queries']:
        assert q['size_query_result']==0 and q['available']
        if q.get('captured'):
            assert q['complete_query_result']==0
            p=w/'query-payloads'/q['file'];assert seal(p)=={'bytes':q['returned_bytes'],'sha256':q['sha256']}
            payloads[q['type']]=p
        else:
            assert q['requested_bytes']>q['allocation_limit_bytes']==8388608
            uncaptured.append(q)
    assert 0 in payloads and 2 in payloads
    code=elf_text(payloads[1].read_bytes()) if 1 in payloads else None
    raw=payloads[2].read_bytes();assert raw[-1]==0
    lines=raw[:-1].decode('utf-8').splitlines();printed=bytearray();instructions=[]
    for number,line in enumerate(lines,1):
        match=re.search(r';\s*((?:[0-9a-f]{8}\s*)+)$',line)
        if not match:continue
        words=match[1].split();instruction=b''.join(struct.pack('<I',int(word,16)) for word in words)
        instructions.append({'line':number,'text_offset':len(printed),'instruction':line.strip().split(';')[0].strip()})
        printed.extend(instruction)
    if code is not None:assert bytes(printed)==code[:len(printed)]
    clocks=[i for i in instructions if 'MSG_RTN_GET_REALTIME' in i['instruction']]
    assert clocks and (site_count is None or len(clocks)==site_count)
    for clock in clocks:
        number=clock['line'];clock['context']=[{'line':i+1,'text':lines[i]} for i in range(max(0,number-13),min(len(lines),number+18))]
    result['stages'][stage]={'module':seal(module),'result':seal(w/'clock-placement.json'),
        'bridge':seal(w/'placement-query-bridge.json'),'owner':seal(w/'run/owner.json'),
        'query_payloads':{p.name:seal(p) for p in payloads.values()},'uncaptured_queries':uncaptured,'statistics':observation['queries'][0]['statistics'],
        'printed_instruction_bytes':len(printed),'ELF_text_bytes':len(code) if code is not None else None,'unprinted_ELF_tail_bytes':len(code)-len(printed) if code is not None else None,
        'all_printed_instruction_bytes_match_ELF_prefix':True if code is not None else None,'clock_sites':clocks,
        'barrier_instruction_count':sum('s_barrier' in i['instruction'] for i in instructions),
        'native_elapsed_seconds':owner['elapsed_seconds'],'native_peak_sampled_rss_bytes':owner['peak_sampled_rss_bytes'],
        'driver_cold_verified':receipt['driver_cold_verified']}
target=WORK/'placement-analysis.json';assert not target.exists();target.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:{f:v[f] for f in ('printed_instruction_bytes','ELF_text_bytes','unprinted_ELF_tail_bytes','barrier_instruction_count','native_elapsed_seconds','native_peak_sampled_rss_bytes')} for k,v in result['stages'].items()}))
