"""Admit completed compiler/link outputs; retain the original reader failure."""
from pathlib import Path
import hashlib,json,re,struct,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;F=W/'fixture';J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())
assert not (F/'build-admission.json').exists()
p=J(F/'prepared.json');built=J(F/'build-result.json');owner=J(F/'build/owner.json')
assert owner['status']=='terminal' and not owner['passed'] and owner['child_exit']==1 and owner['stop_reason'] is None
assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
assert all(x['absent'] for x in owner['observed_births']) and built['passed'] is False
assert len(built['commands'])==3 and all(x['exit']==0 for x in built['commands'])
assert built['new_ISPC_compilations']==built['new_CPP_compilations']==built['new_links']==1
assert 'proof=associate(exe,' in (F/'build/stderr.log').read_text()
for name,value in p['inputs'].items():assert pin(R/name)==value,name
objects=[built['commands'][0]['product'],built['commands'][1]['product'],p['helper']]
for obj in objects:assert pin(R/obj['path'])=={k:obj[k] for k in ('bytes','sha256')}
sys.path.insert(0,str(R/'out/native-cpu-field-compatibility'));import elf_proof
source=elf_proof.source;anchor='assert kind in (2, 4, 42) and offset + 4 <= symbol[\'size\']'
assert source.count(anchor)==1;source=source.replace(anchor,anchor.replace('(2, 4, 42)','(2, 4, 9, 42)'))
anchor="                    original[offset:offset + 4] = final[offset:offset + 4] = b'\\0' * 4"
assert source.count(anchor)==1
source=source.replace(anchor,"""                    # GNU ld may relax a writer-address GOT load into RIP LEA.
                    # Admit only the exact REX.W MOV/LEA opcode transformation;
                    # independently prove the resulting writer target below.
                    if kind == 42 and original[offset-2] == 0x8b and final[offset-2] == 0x8d:
                        assert _target['name'] and offset >= 3
                        assert original[offset-3] == final[offset-3] and original[offset-3] & 0xf8 == 0x48
                        assert original[offset-1] == final[offset-1] and original[offset-1] & 0xc7 == 5
                        original[offset-2] = final[offset-2]
                        masked.append(dict(offset=offset-2, ELF_type=42, bytes=1,
                                           transformation='exact REX.W GOT MOV to RIP LEA; target checked'))
"""+anchor)
namespace={};exec(compile(source,str(R/'out/native-cpu-field-compatibility/elf_proof.py'),'exec'),namespace)
exe=R/p['executable'];proof=namespace['associate'](exe,[dict(o,path=str(R/o['path'])) for o in objects]);assert proof['passed']
# The newly accepted type is the locally declared 4B PC-relative GOT operand.
header=Path('/usr/include/elf.h');assert re.search(r'#define R_X86_64_GOTPCREL\s+9\s',header.read_text())
Elf=namespace['Elf'];linked=Elf(exe);driver=Elf(R/objects[1]['path']);got=[]
def named(elf,name):
 rows=[s for s in elf.symbols if s['name']==name and s['section']];assert len(rows)==1;return rows[0]
for section,location,kind,addend,target in driver.relocations:
 if kind not in (9,42):continue
 assert addend==-4
 containing=[s for s in driver.symbols if s['info']&15==2 and s['section']==section and s['value']<=location<s['value']+s['size']]
 assert len(containing)==1;original=containing[0];final=named(linked,original['name']);offset=location-original['value']
 disp=struct.unpack_from('<i',linked.function_bytes(final),offset)[0]
 relaxed=kind==42 and linked.function_bytes(final)[offset-2]==0x8d and driver.function_bytes(original)[offset-2]==0x8b
 if not target['name'].startswith('writer_') and not relaxed:continue
 if kind==9:assert target['name'].startswith('writer_')
 address=final['value']+offset+4+disp;destination=named(linked,target['name'])
 if relaxed:
  assert address==destination['value'];method='validated RIP LEA'
 else:
  dynamic=[v for v in linked.relocations if v[1]==address];assert len(dynamic)==1 and dynamic[0][2]==8 and dynamic[0][3]==destination['value']
  sections=[s for s in linked.sections if s[3]<=address and address+8<=s[3]+s[5] and s[1]!=8];assert len(sections)==1
  s=sections[0];stored=struct.unpack_from('<Q',linked.data,s[4]+address-s[3])[0];assert stored==destination['value'];method='GOT word plus dynamic RELATIVE addend'
 got.append(dict(caller=original['name'],relocation_offset=offset,ELF_type=kind,reference_address=address,target=target['name'],target_address=destination['value'],method=method))
writer_targets=[x for x in got if x['target'].startswith('writer_')]
assert len(writer_targets)==6 and len({x['target'] for x in writer_targets})==6
proof['writer_pointer_targets']=got
proof['scope']='Every retained function opcode in fixture writer, driver and unchanged AK math objects matches this existing executable outside declared4B PC32/PLT32/GOTPCREL/REX_GOTPCRELX operands and explicitly checked REX.W MOV-to-LEA relaxations. All six writer pointer operands resolve either directly by validated LEA or through one RELATIVE relocation and stored GOT word to the intended writer symbol. Other relaxed local data references are independently bound to their defined symbols. All other opcode bytes compare exactly. Code association only; no numerical or full-renderer qualification.'
(F/'linked-code-admission.json').write_text(json.dumps(proof,indent=2)+'\n')
text=subprocess.check_output(['ldd','-r',str(exe)],text=True,stderr=subprocess.STDOUT);assert 'undefined symbol' not in text and 'not found' not in text
r=dict(passed=True,executable=dict(path=p['executable'],**pin(exe)),objects=objects,linked_code=pin(F/'linked-code-admission.json'),
 runtime_libraries={str(Path(n).resolve()):pin(Path(n).resolve()) for n in re.findall(r'(/[^\s()]+)',text) if Path(n).is_file()},
 original_owner=pin(F/'build/owner.json'),original_owner_passed=False,ELF_type_authority=dict(path=str(header),**pin(header)),
 new_compilations=0,new_links=0,GOTPCREL9_checked=sum(x['ELF_type']==9 for x in got),writer_pointer_targets_checked=len(writer_targets),
 correction='Both compiler commands and linker exited0. Original reader did not support driver GOTPCREL9 function-pointer operands; original failed owner/build result/source are unchanged. Existing bytes now have full opcode and explicit GOT target joins. No compiler/link/workload repeated.')
(F/'build-admission.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(dict(existing_link_admitted=True,original_owner_still_failed=True,joined_functions=[len(x['functions']) for x in proof['objects']],GOT_targets=len(got))))
