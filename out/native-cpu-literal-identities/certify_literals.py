"""Certify only exact 0/+1/-1 identities against canonical device tables."""
from pathlib import Path
import collections,hashlib,json,re,sys
R=Path.cwd();W=Path(__file__).resolve().parent
assert sys.argv[1:] in ([],['--write'])
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def function(s,name):
 m=re.search(r'^static [^\n]+\b'+name+r'\(',s,re.M);assert m,name
 start=s.index('{',m.end());depth=1;end=start+1
 while depth:depth+=(s[end]=='{')-(s[end]=='}');end+=1
 return s[m.start():end]
def output_class(op,a,b,ac,bc):
 if op==0:return {0:0,0x80000000:0,0x3f800000:1,0xbf800000:-1}.get(a)
 if op==7:return -ac if ac is not None else None
 if op==4:
  if ac==0 or bc==0:return 0
  if ac is not None and bc is not None:return ac*bc
 if op in (2,3):
  if ac==0 and bc is not None:return bc if op==2 else -bc
  if bc==0 and ac is not None:return ac
 if op==5 and bc==1:return ac
 return None
def kind(op,ac,bc):
 if op==4 and (ac is not None or bc is not None):return 'multiply'
 if op in (2,3) and (ac==0 or bc==0):return 'add' if op==2 else 'subtract'
 if op==5 and bc==1:return 'divide_one'
 return None
basis=R/'out/native-cpu-direct-compound/ispc-derivation.json';old=json.loads(basis.read_text())
inputs={str(p.relative_to(R)):pin(p) for p in (basis,Path(__file__).resolve())};stages={};total=0;allcounts=collections.Counter()
for slug,item in old['stages'].items():
 p=R/item['derived_path'];assert pin(p)==item['derived'];source=p.read_text()
 table=R/'bin/linux-gcc/native-cpu-slang-bridge/typed-batch'/(slug+'.program.hpp')
 words=list(map(int,re.findall(r'(\d+)u',table.read_text())))
 inputs[str(p.relative_to(R))]=pin(p);inputs[str(table.relative_to(R))]=pin(table)
 known=function(source,'EvaluateKnownLayer_0');headers=list(re.finditer(r'if\(programBase == tableBase \+ (\d+)U && count == (\d+)U && registers == (\d+)U\)',known))
 families={};counts=collections.Counter()
 for hi,h in enumerate(headers):
  offset,n,registers=map(int,h.groups());assert words[offset:offset+2]==[n,registers]
  branch=known[h.end():headers[hi+1].start() if hi+1<len(headers) else len(known)]
  calls=re.findall(r'case (\d+)U:\s*if\(start != (\d+)U \|\| end != (\d+)U\) return false;\s*EvaluateStaticLayer_(\d+)_(\d+)_0\(row, rowBase, context\);',branch)
  assert calls;family,={int(x[3]) for x in calls};assert family not in families
  inst=offset+3+words[offset+2];bounds=inst+5*n;state=[None]*registers;sites=[];i=0
  for layer,start,end,f,fnlayer in calls:
   layer,start,end,f,fnlayer=map(int,(layer,start,end,f,fnlayer));assert f==family and layer==fnlayer
   assert [start,end]==words[bounds+layer:bounds+layer+2] and start==i
   body=function(source,f'EvaluateStaticLayer_{family}_{layer}_0')
   instructions=re.findall(r'EvaluateOneStatic_0<(\d+)U, (\d+)U>\((\d+)U, (\d+)U, (\d+)U, (\d+)U,',body)
   assert len(instructions)==end-start
   for entry in instructions:
    op,reg,dst,a,b,c=map(int,entry);assert reg==registers and dst<reg
    assert [op,dst,a,b,c]==words[inst+5*i:inst+5*i+5]
    arity={0:0,1:0,2:2,3:2,4:2,5:2,6:1,7:1,8:1,9:3,10:3,11:3}[op]
    assert all(x<reg for x in (a,b,c)[:arity])
    ac=state[a] if arity else None;bc=state[b] if arity>=2 else None
    result=output_class(op,a,b,ac,bc);selected=kind(op,ac,bc)
    if selected:
     sites.append({'instruction':i,'layer':layer,'tuple':[op,dst,a,b,c],'a_class':ac,'b_class':bc,'result_class':result,'kind':selected});counts[selected]+=1
    state[dst]=result;i+=1
  assert i==n
  families[str(family)]={'program_offset':offset,'instructions':n,'registers':registers,'layers':len(calls),'site_count':len(sites),'site_sha256':hashlib.sha256(json.dumps(sites,sort_keys=True,separators=(',',':')).encode()).hexdigest()};total+=n
 stages[slug]={'source':{'path':str(p.relative_to(R)),**pin(p)},'canonical_table':{'path':str(table.relative_to(R)),**pin(table)},'families':families,'identity_sites':dict(sorted(counts.items()))};allcounts.update(counts)
assert total==9380
result={'instructions':total,'stages':stages,'identity_sites':dict(sorted(allcounts.items())),'source_inputs':inputs,
 'class_definition':'For every active lane with valid!=0, high/low/tail/error are exactly the numeric constant0,+1or-1 with zero low/tail/error; zero signs and the precise nonzero valid word are retained. Invalid lanes may carry arbitrary words. Unknown remains unknown. No host floating point arithmetic is performed.',
 'rules':['Only literal constructor bit words for0,+1,-1 create initial classes.','Negation flips the class sign.','Multiply by classifiedzero yields zero if valid; multiply two classifiednonzeros yields their sign product if valid.','Add/Subtract propagate only when one classifiedoperand is zero and the other is classified.','Divide propagates a classified numerator only with classified+1denominator.','Every other result becomes unknown; all source reads precede the destination class update, including reused registers.'],
 'scope':'Source-only certificate. All literal tuples and layer/header guards are matched against independently stored canonical device tables. Classes describe all possible valid words, not assumed physical input values; invalid outputs/zero signs/identity precedence still require faithful device implementations and finite whole-word checks. This is metadata specialization, not host-rounded constants, rational precomputation, omitted instruction calls, a new formula or a performance/qualification result.'}
p=W/'literal-certificate.json'
if p.exists():assert json.loads(p.read_text())==result
else:assert sys.argv[1:]==['--write'];p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'instructions_verified_against_tables':total,'identity_sites':dict(sorted(allcounts.items())),'per_stage':{k:v['identity_sites'] for k,v in stages.items()},'compiled':False}))
