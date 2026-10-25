"""Link the already compiled fixture, or execute its predicate-only mode once."""
from pathlib import Path
import hashlib,json,re,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;P=W/'predicate';J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
mode,=sys.argv[1:];p=J(P/'prepared.json');exe=R/p['executable']
if mode=='link':
 assert not (P/'link-result.json').exists() and not exe.parent.exists();exe.parent.mkdir(parents=True)
 for n,v in p['inputs'].items():assert pin(R/n)==v,n
 subprocess.run(p['argv'],check=True)
 text=subprocess.check_output(['ldd','-r',str(exe)],text=True,stderr=subprocess.STDOUT);assert 'undefined symbol' not in text and 'not found' not in text
 sys.path.insert(0,str(R/'out/native-cpu-field-compatibility'));from elf_proof import associate
 proof=associate(exe,[{**x,'path':str(R/x['path'])} for x in p['objects'][1:]]);assert proof['passed']
 r={'passed':True,'new_compilations':0,'new_links':1,'executable':{'path':p['executable'],**pin(exe)},'linked_code':proof,'runtime_libraries':{str(Path(n).resolve()):pin(Path(n).resolve()) for n in re.findall(r'(/[^\s()]+)',text) if Path(n).is_file()}}
 (P/'link-result.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'passed':True,'new_compilations':0,'new_links':1}))
else:
 assert mode=='run' and not (P/'result.json').exists();built=J(P/'link-result.json');assert pin(exe)=={k:built['executable'][k] for k in ('bytes','sha256')}
 out=subprocess.check_output([str(exe),'--predicate-only'],text=True);r=json.loads(out)
 assert r=={'passed':True,'predicate_calls':1764,'compared_tag_words':42336,'input_words_checked':211680},r
 r.update(executable=built['executable'],linked_code=pin(P/'link-result.json'),scope='Shared ISPC predicate agrees on all49scoped tuples/broadcast/mixed/masks/columns, preserving inactive tags and inputs. Arithmetic and full-program/source/FP/scientific qualification are separate.');(P/'result.json').write_text(json.dumps(r,indent=2)+'\n');print(out,end='')
