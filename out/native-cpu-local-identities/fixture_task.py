"""Build once and compare the actual ISPC identity path with the original scalar."""
from pathlib import Path
import hashlib,json,re,struct,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;F=W/'fixture';J=lambda p:json.loads(p.read_text());P=J(F/'prepared.json');exe=R/P['executable']
def pin(p):
 with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())
def verify():
 for name,v in P['inputs'].items():assert pin(R/name)==v,name
verify();mode,=sys.argv[1:]
if mode=='build':
 assert not exe.parent.exists() and not (F/'build-result.json').exists();exe.parent.mkdir(parents=True)
 result=dict(passed=False,commands=[],new_CPP_compilations=0,new_ISPC_compilations=0,new_links=0,new_symbol_renames=0)
 def save():(F/'build-result.json').write_text(json.dumps(result,indent=2)+'\n')
 save()
 for i,argv in enumerate(P['commands']):
  row=dict(argv=argv,started=True);result['commands'].append(row);save();row['exit']=subprocess.run(argv,cwd=R).returncode
  result[['new_symbol_renames','new_ISPC_compilations','new_links'][i]]+=1;save()
  if row['exit']:raise SystemExit(row['exit'])
  p=Path(argv[-1]);row['product']=dict(path=str(p.relative_to(R)),**pin(p));save()
 sys.path.insert(0,str(R/'out/native-cpu-field-compatibility'));import elf_proof
 source=elf_proof.source;anchor='assert kind in (2, 4, 42) and offset + 4 <= symbol[\'size\']';assert source.count(anchor)==1
 source=source.replace(anchor,anchor.replace('(2, 4, 42)','(2, 4, 9, 42)'));ns={};exec(compile(source,str(R/'out/native-cpu-field-compatibility/elf_proof.py'),'exec'),ns)
 Elf=ns['Elf'];old=Elf(R/P['driver']['path']);new=Elf(R/result['commands'][0]['product']['path']);inverse={v:k for k,v in P['mapping'].items()}
 assert len(inverse)==3 and all(not x['section'] for x in old.symbols if x['name'] in P['mapping'])
 assert sorted(old.descriptor(s,{}) for s in old.symbols)==sorted(new.descriptor(s,inverse) for s in new.symbols)
 allocated={}
 for i,s in enumerate(old.sections):
  if not s[2]&2:continue
  matches=[j for j,t in enumerate(new.sections) if new.names[j]==old.names[i]];assert len(matches)==1;j=matches[0];t=new.sections[j]
  assert (s[1],s[2],s[3],s[5],s[8],s[9])==(t[1],t[2],t[3],t[5],t[8],t[9]) and old.bytes(i)==new.bytes(j)
  allocated[old.names[i]]=dict(bytes=s[5],sha256=hashlib.sha256(old.bytes(i)).hexdigest())
 assert set(allocated)=={new.names[i] for i,s in enumerate(new.sections) if s[2]&2}
 relocs=lambda e,m:sorted((e.section_name(s),p,k,a,e.descriptor(t,m)) for s,p,k,a,t in e.relocations)
 assert relocs(old,{})==relocs(new,inverse)
 objects=[result['commands'][1]['product'],result['commands'][0]['product'],P['helper']]
 proof=ns['associate'](exe,[dict(o,path=str(R/o['path'])) for o in objects]);assert proof['passed'];linked=Elf(exe);targets=[]
 def named(name):
  rows=[s for s in linked.symbols if s['name']==name and s['section']];assert len(rows)==1;return rows[0]
 for o in objects[:2]:
  e=Elf(R/o['path'])
  for section,location,kind,addend,target in e.relocations:
   if kind not in (2,4) or not target['name'].startswith(('local_tagged_','sirius_native_tagged_')):continue
   assert addend==-4
   rows=[s for s in e.symbols if s['info']&15==2 and s['section']==section and s['value']<=location<s['value']+s['size']];assert len(rows)==1
   s=rows[0];final=named(s['name']);offset=location-s['value'];disp=struct.unpack_from('<i',linked.function_bytes(final),offset)[0];dest=named(target['name'])
   assert final['value']+offset+4+disp==dest['value'];targets.append(dict(caller=s['name'],offset=offset,target=target['name'],address=dest['value']))
 assert set(P['mapping'].values())<={t['target'] for t in targets} and set(P['mapping'])<={t['target'] for t in targets}
 proof.update(driver_rename=dict(mapping=P['mapping'],all_allocated_sections_identical=allocated,all_symbols_and_relocations_inverse_equal=True),direct_tagged_targets=targets,scope='Every retained ISPC fixture, renamed original driver and unchanged AM math function opcode matches the actual test ELF outside declared4B PC32/PLT32/GOTPCREL/REX_GOTPCRELX operands. Test-only no-relax link; all three renamed driver targets and three fallback targets resolve to intended linked functions. Original weak predicate is unused by argc1. Finite code association, not whole-renderer qualification.')
 (F/'linked-code.json').write_text(json.dumps(proof,indent=2)+'\n')
 raw=subprocess.check_output(['objdump','-d','--insn-width=16',str(R/objects[0]['path'])],text=True)
 mnemonics=set()
 for line in raw.splitlines():
  fields=line.split('\t')
  if len(fields)<3 or not re.match(r'^\s*[0-9a-f]+:$',fields[0]):continue
  assert fields[1].split()[0]!='62';mnemonics.add(fields[2].strip().split()[0])
 assert not re.search(r'%(?:zmm\d+|k[0-7]|mm[0-7]|st)\b',raw)
 assert not any(m.startswith(('vfm','vfnm','f','aes','sha','tile')) or m in ('ldmxcsr','stmxcsr','vldmxcsr','vstmxcsr') for m in mnemonics)
 text=subprocess.check_output(['ldd','-r',str(exe)],text=True,stderr=subprocess.STDOUT);assert 'undefined symbol' not in text and 'not found' not in text
 result.update(executable=dict(path=P['executable'],**pin(exe)),objects=objects,linked_code=pin(F/'linked-code.json'),fixture_mnemonics=sorted(mnemonics),runtime_libraries={str(Path(n).resolve()):pin(Path(n).resolve()) for n in re.findall(r'(/[^\s()]+)',text) if Path(n).is_file()})
 verify();result['passed']=True;save();print(json.dumps(dict(build_passed=True,new_ISPC_compilations=1,new_CPP_compilations=0,new_links=1)))
else:
 assert mode=='run' and not (F/'result.json').exists();built=J(F/'build-result.json');assert built['passed'] and pin(exe)=={k:built['executable'][k] for k in ('bytes','sha256')}
 output=subprocess.check_output([str(exe)],text=True);r=json.loads(output);assert r['passed'] and r['calls']==555660 and r['compared_words']==66679200 and r['tag_words']==13335840 and r['active_lane_results']==1666980
 verify();r.update(executable=built['executable'],linked_code=pin(F/'linked-code.json'),changed_operator_calls=333396,scope='Actual shared ISPC integer identity-store body before original AM arithmetic, called by unchanged original-scalar reference driver. All scratch/tag words, inactive masks and a/b/d aliases, signed zeros/noncanonical valid/nonfinite/refusal, three metadata formats checked. Full9380raw/provider/C8/benefit/frame and scientific/release acceptance are separate.')
 (F/'result.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
