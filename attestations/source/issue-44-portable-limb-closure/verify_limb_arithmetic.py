"""Validate exact mathematical intervals separately from baseline-word equality."""
from collections import Counter
from pathlib import Path
import hashlib
import importlib.util
import json
import resource
import struct
import sys

sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
resource.setrlimit(resource.RLIMIT_FSIZE,(64*1024*1024,64*1024*1024))
root=Path.cwd()
work=Path(__file__).resolve().parent
old=root/'attestations/software-vulkan/41fe5f2/portable-denominator-preparation/files/out/portable-denominator-preparation'
spec=importlib.util.spec_from_file_location('independent_reference',root/'tests/support/portable_binary32/reference.py')
ref=importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
fixture=json.loads((work/'limb-arithmetic-fixture.json').read_text())
assert len(sys.argv)==3
producer,mode=sys.argv[1:]
assert producer in ('baseline','candidate')
assert mode in ('integer','normal','normal_transport','normal_endpoint')
folder=work/('gpu-limb-arithmetic-'+producer+'-'+mode)
data=(folder/'actual.bin').read_bytes()
assert len(data)==fixture['complete_words_each']*4
rows=list(struct.iter_unpack('<8I',data))
assert len(rows)==fixture['records']
wire=(work/'limb-arithmetic-input.bin').read_bytes()
assert len(wire)==4+64*len(rows)
assert hashlib.sha256(wire).hexdigest()==fixture['input']['sha256']
counts=Counter()
opcounts={str(i):Counter() for i in range(4)}

def interval(words):
    values=[ref.exact(v) for v in words]
    assert all(v is not None for v in values) and values[3]>=0
    center=sum(values[:3])
    return center-values[3],center+values[3]

def transform(a,b):
    # Reused original unit-denominator witness: each unchanged FastTwoSum
    # rounds on the independent Fraction lattice, without production PB32.
    if (a&ref.MASK)<(b&ref.MASK): a,b=b,a
    am,bm=a&ref.MASK,b&ref.MASK
    if mode in ('normal_transport','normal_endpoint') and bm==0 and am<ref.INF:
        return (a if am else a&b),b
    high=ref.add(a,b)
    displaced=ref.add(high,a^ref.SIGN)
    low=ref.add(b,displaced^ref.SIGN)
    if mode!='integer' and 0x0c000000<=am<=0x7e800000 and 0x0c000000<=bm<=0x7e800000:
        high=0 if high&ref.MASK==0 else high
        low=0 if low&ref.MASK==0 else low
    return high,low

def unit_pack(a):
    hi,lo,tail,error=a
    x,y=transform(lo,tail)
    z,w=transform(hi,x)
    x,carry=transform(w,y)
    high,w=transform(z,x)
    low,tail=transform(w,carry)
    assert sum(ref.exact(v) for v in (high,low,tail))==sum(ref.exact(v) for v in a[:3])
    return (high,low,tail,error,1,0,0,0)

for index,(case,row) in enumerate(zip(fixture['cases'],rows,strict=True)):
    op=case['op']
    expected_wire=(op,*case['a'],*case['b'],*case['flags'],case['x'],0,0,0,0)
    assert struct.unpack_from('<16I',wire,4+64*index)==expected_wire
    assert row[4] in (0,1) and row[5:]==(0,0,0),(index,'status/padding')
    required=case['required_status']
    if required>=0:
        assert row[4]==required,(index,case['family'],'required admission/refusal',row)
        counts['required_valid' if required else 'required_refusal']+=1
        opcounts[str(op)]['required_valid' if required else 'required_refusal']+=1
    if not row[4]:
        assert row==(0,)*8,(index,'incomplete canonical refusal')
        counts['refused']+=1
        opcounts[str(op)]['refused']+=1
        continue
    bound=0x7d800000 if op==3 else 0x7b800000
    assert all((v&ref.MASK)<=bound for v in row[:4])
    if case['contract']=='raw_compatibility':
        assert op==3 and case['a'][3]&ref.SIGN
        counts['raw_compatibility_valid']+=1
        continue
    assert not row[3]&ref.SIGN
    al,ah=interval(case['a'])
    if op==3:
        assert case['flags'][0]!=0 and (case['x']&ref.MASK)<=0x7b800000
        x=ref.exact(case['x'])
        assert x is not None
        lower,upper=al+x,ah+x
    else:
        bl,bh=interval(case['b'])
        if op==0:
            assert not bl<=0<=bh,(index,'admitted zero-containing denominator')
            corners=[a/b for a in (al,ah) for b in (bl,bh)]
            lower,upper=min(corners),max(corners)
        elif op==1:
            assert all(case['flags'])
            lower,upper=al+bl,ah+bh
        else:
            assert all(case['flags'])
            corners=[a*b for a in (al,ah) for b in (bl,bh)]
            lower,upper=min(corners),max(corners)
    result_lower,result_upper=interval(row[:4])
    assert result_lower<=lower and upper<=result_upper,(
        index,case['family'],'independent exact interval not enclosed',row)
    if op==0 and case['contract']=='zero':
        assert result_lower==result_upper==0
    if op==0 and case['contract']=='self':
        assert result_lower<=1<=result_upper
    if op==0 and case['contract']=='unit':
        assert row==unit_pack(case['a']),(index,'unit denominator complete pack identity')
        counts['unit_complete_pack_identities']+=1
    if case['family']=='connected_full_three_limb_exact_cancellation':
        assert result_lower==result_upper==0,(index,'exact cancellation')
    counts['mathematical_domain_valid']+=1
    opcounts[str(op)]['mathematical_domain_valid']+=1

assert counts['mathematical_domain_valid'] and counts['refused'] and counts['required_valid'] and counts['required_refusal']
for op,c in opcounts.items():
    assert c['mathematical_domain_valid'] and c['refused'] and c['required_valid'] and c['required_refusal'],(op,c)
assert counts['unit_complete_pack_identities']==42
recovery_count=0
for group in fixture['recovery_groups']:
    assert len(group)>1
    values=[rows[i] for i in group]
    assert all(row==values[0] for row in values) and values[0][4]==1
    recovery_count+=len(group)
original_words_equal=False
if mode!='normal_endpoint':
    original_data=(old/('gpu-triple-baseline-'+mode)/'actual.bin').read_bytes()
    assert len(original_data)==fixture['original_division_prefix']*32
    assert data[:len(original_data)]==original_data,'Original5302 division outputs changed'
    original_words_equal=True
if producer=='candidate':
    baseline=work/('gpu-limb-arithmetic-baseline-'+mode)/'actual.bin'
    assert data==baseline.read_bytes(),'Candidate changed a complete original result word'
owner=json.loads((folder/'owner.json').read_text())
assert owner['passed'] and owner['status']=='terminal' and not owner['cleanup_errors']
assert owner['remaining_owned_processes']==[] and owner['whole_inputs_unchanged'] and owner['source_unchanged']
logs=[json.loads(line) for line in (folder/'stdout.log').read_text().splitlines() if line.startswith('{')]
assert any(r.get('driver_id')==13 and r.get('driver')=='llvmpipe' and r.get('RTE32') for r in logs)
done=next(r for r in logs if r.get('completed_dispatches')==1)
assert done['actual_resident_bytes']==fixture['requested_device_bytes']<=8388608
assert done['untouched_complement_words']==done['malformed_wire_words']==0
assert done['observed_words']==fixture['complete_words_each']
assert done['capture_only']==(producer=='baseline')
assert done['mismatches']==(None if producer=='baseline' else 0)
report={'passed':True,'producer':producer,'mode':mode,'records':len(rows),
    'complete_words':fixture['complete_words_each'],'counts':dict(counts),
    'per_operation':{op:dict(c) for op,c in opcounts.items()},
    'recovery_markers_complete_words_equal':recovery_count,
    'all_valid_results_within_mathematical_preconditions_enclose_independent_exact_input_intervals':True,
    'raw_negative_radius_rows_are_compatibility_only':3,
    'original5302division_words_unchanged':original_words_equal,
    'complete_baseline_candidate_words_equal':producer=='candidate',
    'output_bytes':len(data),'output_sha256':hashlib.sha256(data).hexdigest(),
    'scope':'Finite actual accumulator and connected arithmetic/refusal/recovery. '
        'Private2^124 and public2^120 contracts checked separately. '
        'Application joins/integrated stream/science/performance/fullqualification remain pending.'}
output=work/('limb-arithmetic-check-'+producer+'-'+mode+'.json')
assert not output.exists()
output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
