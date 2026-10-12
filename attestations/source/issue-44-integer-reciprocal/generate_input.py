"""Independent full-word fixture plus final quotient boundary probes."""
from pathlib import Path
import hashlib, importlib.util, json, struct, sys
sys.dont_write_bytecode=True
r=Path.cwd();w=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('independent_reference',r/'tests/support/portable_binary32/reference.py')
ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref);ref.selfcheck()
p=r/'bin/linux-gcc/tests/backend/portable_binary32_reference.bin';raw=p.read_bytes()
magic,version,count=struct.unpack_from('<3I',raw)
assert (magic,version,count)==(0x32334250,1,174268) and len(raw)==12+24*count
records=list(struct.iter_unpack('<6I',raw[12:]));seen=set();extra=[]
# Divisors straddle the integer reciprocal's quotient boundaries. Operand
# numerators straddle independently constructed multiples of d/256.
divisors={0x800000,0x800001,0xfffffe,0xffffff}
for reciprocal_quotient in range(256,512):
 d=0xffffffff//reciprocal_quotient
 divisors.update(v for v in (d-1,d,d+1) if 0x800000<=v<0x1000000)
for d in sorted(divisors):
 for k in (128,129,255,256,257,383,384,385,510,511):
  floor_a=k*d//256
  for a in (floor_a-1,floor_a,floor_a+1):
   if not 0x800000<=a<0x1000000:continue
   for index,(ea,eb) in enumerate(((127,127),(1,254),(254,1),(1,2))):
    aw=(ea<<23)|(a-0x800000)|((index&1)<<31)
    bw=(eb<<23)|(d-0x800000)|((index>>1)<<31)
    if (aw,bw) in seen:continue
    seen.add((aw,bw));extra.append((3,aw,bw,ref.scalar_expected(3,aw,bw),0,1))
rounding_subset=0
for d in (0x800000,0x800001,0xfffffe,0xffffff):
 for center in (d,3*d//2):
  for a in (center-1,center,center+1):
   if not 0x800000<=a<0x1000000:continue
   for ea,eb in ((1,151),(1,150),(1,127),(254,127),(254,126),(2,151)):
    for sa in (0,0x80000000):
     for sb in (0,0x80000000):
      aw=(ea<<23)|(a-0x800000)|sa;bw=(eb<<23)|(d-0x800000)|sb
      if (aw,bw) in seen:continue
      seen.add((aw,bw));extra.append((3,aw,bw,ref.scalar_expected(3,aw,bw),0,1));rounding_subset+=1
assert 0<len(extra)<80000
records+=extra;assert len(records)<=262143
inputs=[len(records)];expected=[]
for op,a,b,high,low,valid in records:
 inputs.extend((op,a,b,0));expected.extend((high,low,valid,0))
(w/'input-v2.bin').write_bytes(struct.pack('<%dI'%len(inputs),*inputs))
(w/'expected-v2.bin').write_bytes(struct.pack('<%dI'%len(expected),*expected))
assert (len(inputs)+len(expected))*4<=8*1024*1024
identity=lambda p:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
report={'original_exact_records':count,'extra_division_boundary_records':len(extra),'total_records':len(records),'distinct_normalized_divisors':len(divisors),'boundary_quotient_integers':[128,129,255,256,257,383,384,385,510,511],'exponent_pairs':[[127,127],[1,254],[254,1],[1,2]],'four_sign_combinations_distributed_across_primary_frames':True,'small_rounding_subset_all_signs_per_frame':True,'rounding_subset_records':rounding_subset,'rounding_subset_exponent_pairs':[[1,151],[1,150],[1,127],[254,127],[254,126],[2,151]],'expected_word_method':'existing independent exact Fraction oracle, unchanged','device_allocation_limit_bytes':8*1024*1024,'requested_bytes':(len(inputs)+len(expected))*4,'original_fixture':identity(p),'input':identity(w/'input-v2.bin'),'expected':identity(w/'expected-v2.bin')}
(w/'fixture-v2.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
