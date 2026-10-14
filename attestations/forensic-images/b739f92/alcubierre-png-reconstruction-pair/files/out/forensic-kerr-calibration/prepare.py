from pathlib import Path
import copy,hashlib,json,os,re,shutil,subprocess
root=Path.cwd();packet=root/'out/forensic-reconstruction-b739f92';work=root/'out/forensic-kerr-calibration';assert not work.exists()
work.mkdir()
def ident(p):
 b=p.read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def write(p,o):p.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n')
prep=json.loads((packet/'preparation.json').read_text());oldpins=json.loads((packet/'pins.json').read_text())['files']
for name,expected in oldpins.items():assert ident(Path(name))==expected,name
for name in ('configure-renderers-scaffold','build-renderers-scaffold','build-comparator','comparator-controls'):
 o=json.loads((packet/name/'owner.json').read_text());assert o['passed'] and o['status']=='terminal' and o['child_exit']==0 and not o['remaining_owned_processes'] and not o['cleanup_errors'] and o['whole_inputs_unchanged'] and o['source_unchanged']
build=root/'bin/linux-gcc/forensic-reconstruction-b739f92/scaffold-retry/renderers';case=next(c for c in prep['cases'] if c['name']=='kerr-png')
outputs=root/'renders/forensic-reconstruction-b739f92';output=outputs/'legacy/kerr-png.png';assert not output.exists();output.parent.mkdir(parents=True,exist_ok=True)
for name in ('owned_guard.py','run_owned.py','run_action.py'):shutil.copy2(packet/name,work/name)
shutil.copy2(Path(__file__),work/'prepare.py')
tools=copy.deepcopy(prep['tools']);deps=set();producers={}
for role,filename in [('legacy','Sirius'),('ported','sirius')]:
 p=build/'bin'/filename;assert p.is_file();ldd=subprocess.check_output(['/usr/bin/ldd',str(p)],text=True);assert 'not found' not in ldd
 (work/(role+'-ldd.txt')).write_text(ldd)
 for path in re.findall(r'(?<!\S)(/[^\s()]+)',ldd):
  q=Path(path);assert q.is_file();deps.add(str(q.resolve(strict=True)))
 (work/(role+'-elf-dynamic.txt')).write_bytes(subprocess.check_output(['/usr/bin/readelf','-d',str(p)]))
 producers[role]={'path':str(p),**ident(p),'ldd':str(work/(role+'-ldd.txt'))}
for name in ('/usr/bin/strace','/usr/bin/env','/usr/bin/readelf','/usr/bin/ldd'):
 tools[name]={'resolved_path':str(Path(name).resolve(strict=True)),**ident(Path(name))}
files=[*(packet/name for name in ('pins.json','actions.json','scaffold-overlay.json','preparation.json','source-manifest.json')),
       *(packet/'build-renderers-scaffold'/name for name in ('owner.json','receipt.json')),
       *(packet/'comparator-controls'/name for name in ('owner.json','receipt.json')),
       build/'CMakeCache.txt',build/'compile_commands.json',build/'build.ninja',
       *(Path(p['path']) for p in producers.values()),*(Path(x) for x in deps),
       root/'bin/linux-gcc/forensic-reconstruction-b739f92/comparator/image_compare']
manifest={'evidence_kind':'NEW_RECONSTRUCTION','runner_revision':prep['runner_revision'],'source_export_revision':prep['source_export_revision'],'legacy_runtime_origin':prep['legacy_runtime_origin'],'source_scaffold_overlay':prep['source_scaffold_overlay'],'source_packet':str(packet),'tools':tools,'producers':producers,'dynamic_libraries':{name:ident(Path(name)) for name in sorted(deps)},'case':case,'role':'legacy','output':str(output),'cwd':prep['runtime_cwd'],'budget_scope':'First full-size original-config Kerr PNG under a one-shot300s diagnostic cap to observe actual progress/cost; incomplete at deadline cannot satisfy the reconstructed comparison or original tape identity. No scaled substitute.','limits':{'seconds':300,'rss_kib':4194304,'reserve_kib':2097152},'historical_tape_identity_claimed':False,'scientific_or_release_qualification_claimed':False,'interpretation_limits':['No pair comparison until both complete original-size outputs are decoded','Historical temperature fields do not reach CPU tracer; normalized original radial law retained','Historical observer azimuth is ignored but selected0 matches','Auto worker count must be observed','Asset open and successful6000x3000 load must be verified; grey fallback is failure','strace is diagnostic instrumentation, not a performance benchmark'],'trace_tool':ident(Path('/usr/bin/strace'))}
write(work/'preparation.json',manifest)
argv=['/usr/bin/env','-C',prep['runtime_cwd'],'LC_ALL=C','TZ=UTC','PATH=/usr/bin:/bin','GIT_CEILING_DIRECTORIES='+str(root/'out'),'GIT_CONFIG_NOSYSTEM=1','/usr/bin/strace','-f','-qq','-yy','-e','trace=execve,openat,newfstatat,mmap','-o','{folder}/runtime.trace','--',*case['future_render_argv_review_only']['legacy']]
allfiles=set(list(oldpins)+[str(p) for p in files]+list(tools)+[str(p) for p in work.iterdir() if p.is_file()]+[str(work/'actions.json'),str(work/'pins.json')])
action={'revision':prep['runner_revision'],'kind':'check','expected_exit':0,'limits':manifest['limits'],'scope':manifest['budget_scope'],'argv':argv,'frozen_inputs':sorted(allfiles)}
write(work/'actions.json',{'legacy-kerr-original-size':action})
pins={str(p):ident(p) for p in [Path(x) for x in allfiles if x!=str(work/'pins.json')]}
write(work/'pins.json',{'files':pins,'file_count':len(pins),'file_records_sha256':hashlib.sha256(json.dumps(pins,sort_keys=True,separators=(',',':')).encode()).hexdigest()})
print(json.dumps({'prepared':True,'pins':len(pins),'producers':producers,'libraries':len(deps),'argv':argv,'limits':manifest['limits'],'output':str(output)}))
