from pathlib import Path
import hashlib,json,re,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;p=json.loads((W/'prepared-domain.json').read_text())
def pin(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
def verify():
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==p['revision'] and not subprocess.check_output(['git','status','--porcelain']).strip()
 for n,v in p['inputs'].items():assert pin(R/n)==v,n
verify();assert not (W/'domain-build-result.json').exists()
for args in [p['compile'],p['link']]:subprocess.run(args,check=True)
exe=R/p['executable'];text=subprocess.check_output(['ldd','-r',str(exe)],text=True,stderr=subprocess.STDOUT);assert 'not found' not in text and 'undefined symbol' not in text
sys.path.insert(0,str(R/'out/native-cpu-field-compatibility'));from elf_proof import associate
proof=associate(exe,[{**p['helper'],'path':str(R/p['helper']['path'])}]);assert proof['passed']
r={'passed':True,'new_fixture_CPP_compilations':1,'new_links':1,'new_arithmetic_compilations':0,'helper':p['helper'],'executable':p['executable'],'executable_identity':pin(exe),'linked_math':proof,'runtime_libraries':{str(Path(n).resolve()):pin(Path(n).resolve()) for n in re.findall(r'(/[^\s()]+)',text) if Path(n).is_file()}}
verify();(W/'domain-build-result.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'passed':True,'arithmetic_compilations':0,'fixture_compilations':1,'links':1}))
