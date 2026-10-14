from pathlib import Path
import json,os,hashlib,importlib.util,sys,time,re
work=Path(sys.argv[1]).resolve(strict=True);name=sys.argv[2];case=work/name
owner=json.loads((case/'owner.json').read_text());manifest=json.loads((work/'preparation.json').read_text());action=json.loads((work/'actions.json').read_text())[name];role=action.get('role',manifest.get('role'));assert role in ('legacy','ported');exe_expected=manifest['producers'][role]
spec=importlib.util.spec_from_file_location('guard',work/'owned_guard.py');g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
rows=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  before=g.proc(int(p.name))
  if before['session']!=owner['child_pid']:continue
  exe=os.readlink(p/'exe')
  if exe!=exe_expected['path']:continue
  raw=(p/'maps').read_text();cmd=(p/'cmdline').read_bytes().split(b'\0');cwd=os.readlink(p/'cwd');data=(p/'exe').read_bytes()
  environment={}
  for entry in (p/'environ').read_bytes().split(b'\0'):
   if b'=' not in entry:continue
   key,value=entry.split(b'=',1);key=key.decode()
   if key.startswith(('SIRIUS_','GTEST_','CTEST_','VK_','MESA_','LP_','GALLIVM_','LD_','LC_','OMP_','GIT_')) or key in ('PATH','TZ'):
    environment[key]=value.decode(errors='strict')
  after=g.proc(int(p.name));assert before['start_ticks']==after['start_ticks']
  seen={}
  role_libraries={str(Path(x).resolve(strict=True)) for x in re.findall(r'(?<!\S)(/[^\s()]+)',Path(exe_expected['ldd']).read_text())}
  assert role_libraries <= set(manifest['dynamic_libraries'])
  for path in sorted(role_libraries):
   pin=manifest['dynamic_libraries'][path]
   mappings=[line for line in raw.splitlines() if line.endswith(' '+path)];assert mappings,path
   b=Path(path).read_bytes();assert len(b)==pin['bytes'] and hashlib.sha256(b).hexdigest()==pin['sha256'];st=Path(path).stat();seen[path]={'file':pin,'mappings':mappings,'stat_inode':st.st_ino,'stat_dev':st.st_dev}
  row={'before':before,'after':after,'utc_ns':time.time_ns(),'exe':exe,'exe_bytes':len(data),'exe_sha256':hashlib.sha256(data).hexdigest(),'cwd':cwd,'argv':[x.decode() for x in cmd if x],'libraries':seen,'declared_role_library_count':len(role_libraries),'other_role_input_libraries_not_required_to_be_mapped':sorted(set(manifest['dynamic_libraries'])-role_libraries),'observed_relevant_environment':environment,'environment_scope':'Only documented renderer/config/provider/loader/test/locale/thread/Git selectors; unrelated environment and credentials deliberately excluded','scope':'One stable-birth live process snapshot; map paths/whole backing-file hashes observed, no general provider admission claim'}
  assert row['exe_sha256']==exe_expected['sha256'];rows.append(row);(case/'live-maps.txt').write_text(raw)
 except (PermissionError,FileNotFoundError,ProcessLookupError):continue
assert len(rows)==1,('expected one live renderer',len(rows),owner)
p=case/'live-bindings.json';assert not p.exists();p.write_text(json.dumps(rows,indent=2)+'\n');print({'role':role,'pid':rows[0]['before']['pid'],'libraries':len(rows[0]['libraries']),'scoped_environment_captured':True})
