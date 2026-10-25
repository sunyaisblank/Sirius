"""Admit the completed link with the intended objects; retain the failed owner."""
from pathlib import Path
import hashlib,json,re,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;P=W/'predicate';J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
assert not (P/'link-result.json').exists();p=J(P/'prepared.json');o=J(P/'link/owner.json')
assert o['status']=='terminal' and not o['passed'] and o['child_exit']==1 and o['stop_reason'] is None
assert o['source_unchanged'] and o['whole_inputs_unchanged'] and not o['remaining_owned_processes'] and not o['cleanup_errors']
assert o['observed_births'] and all(x['absent'] for x in o['observed_births'])
assert o['wall_seconds']<=300 and o['peak_sampled_rss_kib']<=4194304 and o['minimum_sampled_available_kib']>=2097152
stderr=(P/'link/stderr.log').read_text();assert 'proof=associate(exe,' in stderr and stderr.rstrip().endswith('AssertionError')
assert not (P/'link/stdout.log').read_text().strip()
for name,v in p['inputs'].items():assert pin(R/name)==v,name
objects=[x for x in p['objects'] if x['path'].endswith(('/retained8.o','/predicate.o'))];assert len(objects)==2
sys.path.insert(0,str(R/'out/native-cpu-field-compatibility'));from elf_proof import associate
exe=R/p['executable'];proof=associate(exe,[{**x,'path':str(R/x['path'])} for x in objects]);assert proof['passed']
predicate=next(R/x['path'] for x in objects if x['path'].endswith('/predicate.o'))
raw=subprocess.check_output(['objdump','-d','--insn-width=16',str(predicate)],text=True);mnemonics=set()
for line in raw.splitlines():
 fields=line.split('\t')
 if len(fields)>=3 and re.match(r'^\s*[0-9a-f]+:$',fields[0]):
  assert fields[1].split()[0]!='62';mnemonics.add(fields[2].strip().split()[0])
forbidden=[x for x in mnemonics if x.startswith(('vfm','vfnm')) or x in ('ldmxcsr','stmxcsr','vldmxcsr','vstmxcsr') or re.match(r'^v?(?:add|sub|mul|div|sqrt|round|cvt|cmp|ucomi|comi)[a-z]*(?:ps|ss|pd|sd)$',x)]
assert not forbidden and not re.search(r'%(?:zmm\d+|k[0-7]|mm[0-7]|st)\b',raw),forbidden
text=subprocess.check_output(['ldd','-r',str(exe)],text=True,stderr=subprocess.STDOUT);assert 'undefined symbol' not in text and 'not found' not in text
r={'passed':True,'new_compilations':0,'new_links':1,'executable':{'path':p['executable'],**pin(exe)},'linked_code':proof,'runtime_libraries':{str(Path(n).resolve()):pin(Path(n).resolve()) for n in re.findall(r'(/[^\s()]+)',text) if Path(n).is_file()},'original_link_owner':pin(P/'link/owner.json'),'original_link_owner_passed':False,'verification_reader':pin(Path(__file__).resolve()),'predicate_integer_ISA':{'object':pin(predicate),'mnemonics':sorted(mnemonics),'forbidden':forbidden},'correction':'The original link task completed its checked linker invocation and then selected objects[1:] for code association, including driver.o instead of retained8.o. Current reader explicitly selects the unchanged math and predicate objects and joins every retained function. Original failed owner, stderr, preparation and source remain unchanged; no compiler, linker or workload repeated. This result admits the existing link/code association only, not the failed original owner or numerical/scientific qualification.'}
(P/'link-result.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'existing_owned_link_admitted':True,'original_owner_still_failed':True,'repeated_links':0,'joined_functions':[len(x['functions']) for x in proof['objects']]}))
