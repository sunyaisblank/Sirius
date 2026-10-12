"""Independent exact dyadic words and original two-round seed expectations."""
from pathlib import Path
from fractions import Fraction
import hashlib,importlib.util,json,random,struct,sys
sys.dont_write_bytecode=True
r=Path.cwd();w=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('independent_reference',r/'tests/support/portable_binary32/reference.py');ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref);ref.selfcheck()
p=r/'bin/linux-gcc/tests/backend/portable_binary32_reference.bin';raw=p.read_bytes();assert struct.unpack_from('<3I',raw)==(0x32334250,1,174268) and len(raw)==12+24*174268
original=list(struct.iter_unpack('<6I',raw[12:]));edges=[row[2] for row in original[:74]];assert len(set(edges))==74
power=lambda e:Fraction(1<<e) if e>=0 else Fraction(1,1<<-e)
def scale(a,e):
 c=ref.classify(a)
 return ref.NAN if c==4 else a if c in (0,3) else ref.round_exact(ref.exact(a)*power(e),a&ref.SIGN)
def divide_seed(a,b):
 if ref.classify(a)==0:return 0
 aa=abs(ref.exact(a));bb=abs(ref.exact(b));assert bb>0
 ea=ref.exponent_of(aa);eb=ref.exponent_of(bb);e=ea-eb
 if e < -151 or e > 120:return 0
 seed=ref.round_exact((aa/power(ea))/(bb/power(eb)))
 if e < -126:seed=scale(seed,-126);e+=126
 return scale(seed,e)^((a^b)&ref.SIGN)
def sqrt_seed(a):
 aa=ref.exact(a);assert aa>0;e=ref.exponent_of(aa);odd=e%2
 seed=ref.round_sqrt((aa/power(e))*power(odd))
 return scale(seed,(e-odd)//2)
scales=[]
for e in range(-149,128):
 for a in edges:scales.append((14,a,e&0xffffffff,scale(a,e),0,1))
rng=random.Random(0x5349524955530073)
for _ in range(8192):
 a=rng.getrandbits(32);e=rng.randrange(-149,128);scales.append((14,a,e&0xffffffff,scale(a,e),0,1))
seeds=[]
for _ in range(8192):
 a=rng.randrange(0x7b800001)|rng.choice((0,ref.SIGN));b=rng.randrange(1,0x7b800001)|rng.choice((0,ref.SIGN));seeds.append((15,a,b,divide_seed(a,b),0,1))
for _ in range(4096):
 a=rng.randrange(1,0x7b800001);seeds.append((16,a,0,sqrt_seed(a),0,1))
for sa in (0,ref.SIGN):
 for sb in (0,ref.SIGN):
  a=0x00800000|sa;b=0x40fffff6|sb;expected=0x00100000|(sa^sb);assert divide_seed(a,b)==expected;seeds.append((15,a,b,expected,0,1))
assert scale(0x3f000005,-128)==0x00100001 and scale(scale(0x3f000005,-126),-2)==0x00100000
records=original+scales
identity=lambda p:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
report={'original_records':len(original),'scale_edge_cross_records':277*74,'scale_random_records':8192,'scale_exponents':[-149,127],'cpu_records':len(records),'seed_division_records':8192+4,'seed_sqrt_records':4096,'two_rounding_counterexample':{'a':'0x00800000','b':'0x40fffff6','seed':'0x3f000005','separate_rounds':'0x00100000','merged_round':'0x00100001'},'expectation_method':'unchanged independent Fraction lattice/round-sqrt oracle; seed decomposition retains separately required rounding steps, no native arithmetic','source_fixture':identity(p)}
for mode,rows in [('cpu',records),('gpu',records+seeds)]:
 assert len(rows)<=262143;inputs=[len(rows)];expected=[]
 for op,a,b,h,l,v in rows:inputs.extend((op,a,b,0));expected.extend((h,l,v,0))
 assert 4*(len(inputs)+len(expected))<=8*1024*1024
 for name,words in [('input',inputs),('expected',expected)]:
  target=w/(mode+'-'+name+'.bin');target.write_bytes(struct.pack('<%dI'%len(words),*words));report[mode+'_'+name]=identity(target)
 report[mode+'_requested_bytes']=4*(len(inputs)+len(expected));report[mode+'_records']=len(rows)
(w/'fixture.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
