"""Keep the original division corpus and add exact interval/carry partitions."""
from collections import Counter
from pathlib import Path
import hashlib
import importlib.util
import itertools
import json
import random
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
ref.selfcheck()
original=json.loads((old/'triple-fixture.json').read_text())
assert original['records']==5302
cases=[dict(c,op=0,flags=[0,0],x=0) for c in original['cases']]
counts=Counter(c['family'] for c in cases)
markers=[list(original['recovery_indices'])]
PUBLIC=0x7b800000
PRIVATE=0x7d800000
unit=(0x3f800000,0,0,0)
zero=(0,0,0,0)

def power(e):
    assert -149<=e<=124
    return (e+127)<<23 if e>=-126 else 1<<(e+149)

def append(op,a,b,family,status=1,contract='interval',flags=(1,1),x=0):
    assert len(a)==len(b)==4
    assert all(0<=v<=0xffffffff for v in (*a,*b,*flags,x))
    cases.append({'op':op,'a':list(a),'b':list(b),'family':family,
        'required_status':status,'contract':contract,'flags':list(flags),'x':x})
    counts[family]+=1

# Direct fourth-limb loss: terms are separated by 24 bits, so each addition
# has an exact rational reference even when a residual reaches the error radius.
for e in (-77,-30,0,30,60,96,120):
    for signs in itertools.product((0,ref.SIGN),repeat=4):
        terms=[power(e-24*i)|sign for i,sign in enumerate(signs)]
        for missing in range(4):
            selected=[v for i,v in enumerate(terms) if i!=missing]
            append(3,(*selected,0),zero,'direct_three_carries_and_fourth_residual',x=terms[missing])

# Zero signs and adjacent subnormal/normal/fast-sum/public/private boundaries.
for signs in itertools.product((0,ref.SIGN),repeat=4):
    append(3,(*signs[:3],0),zero,'direct_signed_zero_limbs',x=signs[3])
for boundary in (1,0x00800000,0x0c000000,0x3f800000,0x5d800000,PUBLIC,PRIVATE):
    for offset in (-1,0,1):
        magnitude=boundary+offset
        for sign in (0,ref.SIGN):
            word=magnitude|sign
            # Public x bound refuses either sign above 2^120; private input
            # bounds permit larger limbs when the result stays at <=2^124.
            append(3,zero,zero,'direct_public_x_adjacent_boundaries',
                   status=1 if magnitude<=PUBLIC else 0,x=word)
            if magnitude<=PRIVATE:
                append(3,(word,0,0,0),zero,'direct_private_limb_adjacent_boundaries',x=0)
for e in (-149,-148,-127,-126,-104,-103,-102,-25,-24,-1,0,60,120):
    for sign in (0,ref.SIGN):
        v=power(e)|sign
        append(3,(v,v^ref.SIGN,1,1),zero,'direct_large_term_cancellation',x=1)
        append(3,(v,1,1,1),zero,'direct_subnormal_residuals',x=power(max(-149,e-25)))
for e in (-103,-30,0,30,96):
    for parity in (0,1):
        for sign in (0,ref.SIGN):
            append(3,((power(e)+parity)|sign,power(e-24),power(max(-149,e-48)),0),zero,
                   'direct_odd_even_half_ulp_carries',x=power(max(-149,e-72)))
append(3,(0x3f800000,0x3f800000,0x3f800000,0),zero,
       'direct_overlapping_limbs_fourth_residual',x=power(-25))
for negative_radius in (ref.SIGN,ref.SIGN|1,0xbf800000):
    append(3,(*unit[:3],negative_radius),zero,'raw_negative_radius_compatibility',
           -1,'raw_compatibility',x=0)
for flag in (0,1,2,0xffffffff):
    append(3,(0x3f800000,power(-24),power(-48),0),zero,'direct_valid_flag_predicate',
           status=int(flag!=0),flags=(flag,1),x=power(-72))
for position in range(4):
    for bad in (PRIVATE+1,ref.INF,ref.NAN):
        a=list(unit);a[position]=bad
        append(3,a,zero,'direct_nonfinite_or_private_out_of_bound',0,x=0)
for radius in (0,1,power(120),PRIVATE):
    append(3,(*zero[:3],radius),zero,'direct_private_error_bounds',x=0)
append(3,(PRIVATE,0,0,0),zero,'direct_private_high_overflow',0,x=PUBLIC)
append(3,(PRIVATE|ref.SIGN,0,0,0),zero,'direct_private_negative_high_overflow',0,x=PUBLIC|ref.SIGN)

# Addition: each safe input term <=2^116 leaves ample public headroom.
# Product: each limb/radius <=2^50 leaves headroom for all nine terms and
# interval uncertainties. Hard refusal checks below separately retain limits.
for op in (1,2):
    for signs in itertools.product((0,ref.SIGN),repeat=3):
        append(op,(*signs,0),unit,'connected_signed_zero_and_identity')
        append(op,unit,(*signs,0),'connected_identity_and_signed_zero')
    for e in (-149,-127,-126,-103,-60,-30,-1,0,1,30,50):
        for sign in (0,ref.SIGN):
            a=(power(e)|sign,power(max(-149,e-25)),power(max(-149,e-50)),1)
            for b in (unit,(0xbf800000,0,0,0),(0x3fc00000,power(-25),power(-50),1),a):
                append(op,a,b,'connected_scale_three_limb_uncertainty')
    for bad in (PUBLIC+1,ref.INF,ref.NAN):
        for position in range(4):
            a=list(unit);a[position]=bad
            append(op,a,(0x40000000,0,0,0),'connected_invalid_pack',0)
            append(op,(0x40000000,0,0,0),a,'connected_invalid_pack',0)
    for negative_radius in (ref.SIGN,ref.SIGN|1):
        append(op,(*unit[:3],negative_radius),unit,'connected_negative_radius_refusal',0)
    for flags in ((0,1),(1,0),(0,0)):
        append(op,unit,unit,'connected_invalid_flag_refusal',0,flags=flags)

for sign in (0,ref.SIGN):
    append(1,(PUBLIC|sign,0,0,0),(PUBLIC|sign,0,0,0),'connected_public_sum_overflow',0)
    # Magnitude at 2^60 is accepted by product transforms; a next-word value
    # is refused when neither operand triggers the algebraic identity paths.
    append(2,(power(60)|sign,0,0,0),(0x40000000,0,0,0),'connected_product_factor_at_bound')
    append(2,((power(60)+1)|sign,0,0,0),(0x40000000,0,0,0),
           'connected_product_factor_above_bound',0)
for e in (-103,0,60,116):
    append(1,(power(e),power(max(-149,e-24)),power(max(-149,e-48)),0),
           (power(e)|ref.SIGN,power(max(-149,e-24))|ref.SIGN,power(max(-149,e-48))|ref.SIGN,0),
           'connected_full_three_limb_exact_cancellation')

rng=random.Random(0x5349524955530313)
def random_triple(ceiling):
    e=rng.randrange(-120,ceiling+1)
    hi=power(e)|rng.randrange(1<<23)|rng.choice((0,ref.SIGN))
    return (hi,power(max(-149,e-25))|rng.choice((0,ref.SIGN)),
            power(max(-149,e-50))|rng.choice((0,ref.SIGN)),
            rng.choice((0,1,power(max(-149,e-60)))))
for _ in range(2048):
    a=random_triple(116)
    append(3,a,zero,'direct_fixed_seed_carry_cases',x=random_triple(116)[0])
    append(1,a,random_triple(116),'connected_fixed_seed_addition')
    append(2,random_triple(50),random_triple(50),'connected_fixed_seed_product')

# Refusal adjacency and workgroup-crossing repeated inputs catch accidental
# state leaking across operations. Preserve the 5302 original division rows.
added=cases[5302:]
cases=cases[:5302]
for op in (1,2,3):
    group=[]
    for index,c in enumerate(c for c in added if c['op']==op):
        if index%31==0 or c['required_status']==0:
            group.append(len(cases))
            append(op,(0x3fc00000,power(-25),power(-50),1),
                   (0x3fc00000,power(-25),power(-50),1),'added_recovery_marker',
                   contract='recovery',x=power(-75))
        cases.append(c)
        if c['required_status']==0:
            group.append(len(cases))
            append(op,(0x3fc00000,power(-25),power(-50),1),
                   (0x3fc00000,power(-25),power(-50),1),'added_recovery_marker',
                   contract='recovery',x=power(-75))
    markers.append(group)

inputs=[len(cases)]
for case in cases:
    inputs.extend((case['op'],*case['a'],*case['b'],*case['flags'],case['x'],0,0,0,0))
assert len(inputs)==1+16*len(cases)
wire=struct.pack('<%dI'%len(inputs),*inputs)
oldwire=(old/'triple-input.bin').read_bytes()
assert wire[4:4+5302*64]==oldwire[4:]
words=8*len(cases)
requested=4*(len(inputs)+words)
assert requested<=8*1024*1024 and len(cases)<=87381
# The contracts are derived from exact sums and interval arithmetic, not the
# candidate accumulation sequence. Decode only finite valid-domain operands.
for c in cases[5302:]:
    if c['required_status']!=1: continue
    vals=[ref.exact(v) for v in c['a']]
    assert all(v is not None for v in vals) and vals[3]>=0
    if c['op']==3:
        assert ref.exact(c['x']) is not None
    else:
        assert all(ref.exact(v) is not None for v in c['b']) and ref.exact(c['b'][3])>=0

for name,data in (('limb-arithmetic-input.bin',wire),('limb-arithmetic-capture-shape.bin',bytes(words*4))):
    p=work/name
    assert not p.exists()
    p.write_bytes(data)
identity=lambda p:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
report={'records':len(cases),'original_division_prefix':5302,'complete_words_each':words,
    'requested_device_bytes':requested,'partitions':dict(Counter(c['family'] for c in cases)),
    'operations':dict(Counter(str(c['op']) for c in cases)),
    'input':identity(work/'limb-arithmetic-input.bin'),'original_input':identity(old/'triple-input.bin'),
    'original_fixture':identity(old/'triple-fixture.json'),'recovery_groups':markers,'cases':cases,
    'expectation_method':'Independent exact Fraction sum and interval sum/product/quotient enclosure. '
        'Original/candidate complete words are a separate compatibility target. All added safe rows require admission; '
        'nonfinite/known bound/zero flags require canonical refusal. Direct private output uses 2^124 bounds.',
    'scope':'Finite actual accumulator/connected add,multiply,division/refusal/recovery witnesses; '
        'application joins, integrated scientific/stream, performance and full qualification remain separate.'}
(work/'limb-arithmetic-fixture.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ('cases','recovery_groups')}))
