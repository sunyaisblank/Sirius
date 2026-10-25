"""Summarize, preserve, and retire this single terminal diagnostic's scratch."""
from pathlib import Path
from collections import Counter
import hashlib,json,os,subprocess,sys,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent
J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())
def write(p,v):
 with p.open('x') as f:json.dump(v,f,indent=2);f.write('\n')
mode,=sys.argv[1:]
plan=J(W/'plan.json');F=W/'observe';B=R/'attestations/qualification'/plan['source'][:8]/'native-cpu-steady-profile'
if mode=='preserve':
 owner=J(F/'owner.json');result=J(W/'result.json');raw=J(F/'all-threads-observation.json');decoded=J(F/'all-threads-symbols.json')
 assert owner['status']=='terminal' and owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
 assert not owner['remaining_owned_processes'] and not owner['cleanup_errors'] and all(x['absent'] for x in owner['observed_births'])
 assert result['status']=='complete' and result['intentional_stop'] and result['original_command_exit']!=0 and not result['frame_qualification']
 assert result['trigger_progress'][0]>=128 and result['trigger_progress'][1]==16384
 assert raw['status']=='complete' and raw['completed_requested_window']
 assert not raw['errors'] and not raw['incomplete_reasons'] and not raw['cleanup_errors']
 assert decoded['status']=='decoded' and not decoded['errors']
 libraries=Counter();objects=Counter();native=0
 for row in decoded['aggregate']['functions']:
  libraries[row.get('file') or 'unresolved']+=row['samples']
  assert len(row['native_objects'])<=1
  for name in row['native_objects']:objects[name]+=row['samples'];native+=row['samples']
 assert sum(libraries.values())==raw['sample_count']==decoded['aggregate']['samples']
 analysis=dict(source=plan['source'],diagnostic_complete=True,frame_accepted=False,scientific_qualification=False,
  trigger_progress=result['trigger_progress'],capture_samples=raw['sample_count'],stable_observed_threads=len(raw['threads']),
  native_joined_samples=native,samples_by_file=dict(libraries),samples_by_native_object=dict(objects),
  top_functions=decoded['aggregate']['functions'][:20],observer_wall_seconds=result['observer_wall_seconds'],
  observer_peak_rss_kib=result['observer_peak_rss_kib'],owned_wall_seconds=owner['wall_seconds'],
  owned_peak_sampled_rss_kib=owner['peak_sampled_rss_kib'],minimum_sampled_available_kib=owner['minimum_sampled_available_kib'],
  conclusion='This steady body sample places the majority of observed USER CPU PCs in verified native retained arithmetic/execution. Tagged multiply/divide are prominent. It does not support treating inclusive API wait/dispatch windows as dominant isolated CPU cost. Further work should reduce faithfully executed retained arithmetic work or its recurring cost with an explicit equivalence argument; another scheduling/ABI-only or arithmetic-layout variation is not justified by this capture.',
  limits='One3s100Hz USER-only PC window after128 completed tiles,100 stable observed threads. No callchains, kernel/off-CPU attribution, whole-frame causal shares, timing comparison, parameter sweep or completed image. Deliberately stopped at192/16384. Seven samples have unresolved maps; source/whole frozen inputs unchanged, terminal cleanup complete. Prior initialization capture is separate and not pooled.')
 write(W/'analysis.json',analysis)
 paths=sorted(p for p in W.rglob('*') if p.is_file())
 paths+=[R/'out/native-cpu-runtime-observe/all_threads.py',R/'out/native-cpu-runtime-observe/decode_all_threads.py',R/'out/native-cpu-normalization-tags/provider/linked-native-proof.json']
 B.mkdir(parents=True,exist_ok=False);archive=B/'source-and-observation.zip';members={}
 with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
  for p in paths:
   name=str(p.relative_to(R));assert not p.is_symlink()
   members[name]={**pin(p),**{k:getattr(p.stat(),k) for k in ('st_mode','st_mtime_ns')}};z.write(p,name)
 dependencies=[]
 for p in (R/'attestations/qualification/0ce77927/native-cpu-normalization-tags-frame-timeout/actual-products.zip',R/'attestations/qualification/f878a18f/native-cpu-normalized-arithmetic-frame-timeout/actual-products.zip'):
  assert p.is_file();dependencies.append(dict(path=str(p.relative_to(R)),**pin(p)))
 manifest=dict(source=plan['source'],scope=plan['scope'],frame_accepted=False,archive=dict(path=str(archive.relative_to(R)),**pin(archive)),
  members=members,current_runtime_product_recovery=dependencies,analysis=pin(W/'analysis.json'),terminal_owner=pin(F/'owner.json'),
  full_scientific_and_release_qualification=False)
 write(B/'MANIFEST.json',manifest)
 # Independently reopen the finished ZIP and compare every payload to source.
 with zipfile.ZipFile(archive) as z:
  assert z.testzip() is None and set(z.namelist())==set(members)
  for name,v in members.items():
   data=z.read(name);p=R/name
   assert len(data)==v['bytes'] and hashlib.sha256(data).hexdigest()==v['sha256'] and data==p.read_bytes()
   assert z.getinfo(name).external_attr>>16==v['st_mode']
 print(json.dumps(dict(preserved_members=len(members),archive=pin(archive),passed_diagnostic_only=True)))
elif mode=='retire':
 manifest=J(B/'MANIFEST.json');archive=R/manifest['archive']['path']
 assert pin(archive)=={k:manifest['archive'][k] for k in ('bytes','sha256')}
 for p in (B/'MANIFEST.json',archive):
  assert subprocess.check_output(['git','show','@{upstream}:'+str(p.relative_to(R))])==p.read_bytes()
 restored=J(W/'resource-restore.json');cache=R/restored['cache'];dest=R/restored['destination']
 expected=restored['files'];assert {str(p.relative_to(dest)) for p in dest.rglob('*') if p.is_file()}==set(expected)
 for name,value in expected.items():assert os.path.samefile(cache/name,dest/name) and pin(dest/name)==value and pin(cache/name)==value
 temporary=Path(J(W/'actions.json')['observe']['build_environment']['TMPDIR']).parent
 assert temporary.is_relative_to(R/'renders') and not any(p.is_file() or p.is_symlink() for p in temporary.rglob('*'))
 retired=[];allocated=0
 with zipfile.ZipFile(archive) as z:
  assert z.testzip() is None
  for p in sorted(F.iterdir()):
   assert p.is_file() and not p.is_symlink();name=str(p.relative_to(R));value=manifest['members'][name]
   assert pin(p)=={k:value[k] for k in ('bytes','sha256')} and z.read(name)==p.read_bytes()
   assert not subprocess.check_output(['git','ls-files','--',name]).strip()
   allocated+=p.stat().st_blocks*512;retired.append(dict(path=name,**pin(p)));p.unlink()
 F.rmdir()
 for name in expected:(dest/name).unlink()
 for base in (dest,temporary):
  for p in sorted(base.rglob('*'),key=lambda p:len(p.parts),reverse=True):p.rmdir()
  base.rmdir()
 write(W/'retirement.json',dict(source=plan['source'],verified_pushed_recovery=manifest['archive'],
  retired_diagnostic_files=retired,gross_diagnostic_allocated_bytes=allocated,resource_aliases_removed=len(expected),
  resource_payload_bytes_removed=0,temporary_render_tree_removed=str(temporary.relative_to(R)),
  remaining_owned_processes=[],net_drive_gain_claimed=False))
 print(json.dumps(dict(expanded_diagnostic_files_retired=len(retired),gross_allocated_bytes=allocated,resource_aliases_removed=32,unique_resources_preserved=True)))
else:raise SystemExit('preserve or retire')
