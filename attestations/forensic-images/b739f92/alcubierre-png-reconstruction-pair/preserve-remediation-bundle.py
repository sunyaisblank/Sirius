from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,tempfile
root=Path.cwd();kind=sys.argv[1]
def identity(p):
 if p.is_symlink():
  s=os.readlink(p);b=s.encode();return {'type':'symlink','target':s,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
 b=p.read_bytes();return {'type':'file','bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'mode':p.stat().st_mode&0o777}
def write(p,x):p.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
source=subprocess.check_output(['git','rev-parse','HEAD']).decode().strip();assert source=='f5d649bf5067d1d68dc9215462f37830122f0540' and not subprocess.check_output(['git','status','--porcelain']).strip()
if kind=='translation':
 w=root/'out/macos-moltenvk-unforced';bundle=root/'attestations/native-vulkan/4bcd987/macos-unforced-translator-match';tag='evidence/issue-48-unforced-translator-match'
 accepted=json.loads((w/'observations/translation-byte-acceptance.json').read_text());assert accepted['passed'] and accepted['stock']['sha256']=='b2d7daec9928356de4ca24a724d135adc59685a330b043e4fa55f2fe2f9b5550' and accepted['candidate']['sha256']=='80b33059934cc77a685c5738b9413c94379e45f1b18d0e6635a6432ff0407921'
 summary={'schema':'sirius-unforced-translator-result-v1','runner_revision':source,'product_source':'4bcd9875cd9308deabed06a5fa2f17d2321b0cbc','source_freeze':identity(w/'input-freeze.json'),'byte_acceptance':accepted,'content_id':json.loads((w/'content-id.json').read_text())['content_id'],'native_provider_built_or_executed':False,'qualification_claimed':False,'owners':{},'lineage':json.loads((w/'candidate-translation/lineage.json').read_text())}
 for name in ('stock-translation','candidate-translation'):
  o=json.loads((w/name/'owner.json').read_text());assert o['passed'] and o['status']=='terminal' and o['child_exit']==0 and o['source_unchanged'] and o['whole_inputs_unchanged'] and not o['remaining_owned_processes'] and not o['cleanup_errors'];assert all(x['absent'] for x in o['observed_births'])
  summary['owners'][name]={k:o[k] for k in ('wall_seconds','peak_sampled_rss_kib','stop_reason','observed_births','limits')}
 mapping=[('out/macos-moltenvk-unforced','files/out/macos-moltenvk-unforced'),('bin/linux-gcc/macos-moltenvk-unforced','files/bin/linux-gcc/macos-moltenvk-unforced')]
 message='Preserve exact stock and annotation-only translator joins (#48)\n'
elif kind in ('forensic-png','forensic-exr','forensic-alc'):
 w=root/('out/forensic-alcubierre-completion' if kind=='forensic-alc' else 'out/forensic-kerr-exr' if kind=='forensic-exr' else 'out/forensic-kerr-calibration');label='alcubierre-png' if kind=='forensic-alc' else 'kerr-exr' if kind=='forensic-exr' else 'kerr-png';bundle=root/('attestations/forensic-images/b739f92/'+label+'-reconstruction-pair');tag='evidence/issue-100-'+label+'-reconstruction-pair'
 summary=json.loads((w/'pair-result.json').read_text());assert summary['reconstruction_pair_complete'] and summary['comparison']['equal'] and summary['comparison']['mismatched_samples']==0
 mapping=[('out/forensic-reconstruction-b739f92','files/out/forensic-reconstruction-b739f92'),('out/forensic-kerr-calibration','files/out/forensic-kerr-calibration'),('out/forensic-kerr-ported','files/out/forensic-kerr-ported'),('bin/linux-gcc/forensic-reconstruction-b739f92','files/bin/linux-gcc/forensic-reconstruction-b739f92'),('renders/forensic-reconstruction-b739f92','files/renders/forensic-reconstruction-b739f92')]
 if kind in ('forensic-exr','forensic-alc'):mapping.append(('out/forensic-kerr-exr','files/out/forensic-kerr-exr'))
 if kind=='forensic-alc':
  assert summary['comparison']['width']==1920 and summary['comparison']['height']==1080 and summary['comparison']['compared_samples']==8294400 and summary['rays_per_renderer']==132710400
  assert not summary['historical_tape_identity_claimed'] and not summary['current_physics_or_scientific_or_release_qualification_claimed']
  for name in ('legacy-alcubierre-original-size','ported-alcubierre-original-size'):
   r=json.loads((w/name/'result.json').read_text());o=json.loads((w/name/'owner.json').read_text())
   assert r['classification']=='COMPLETE_ORIGINAL_SIZE_SINGLE_RENDER' and r['complete_original_size_render'] and r['self_decode']['passed'] and not r['proof_gaps']
   assert o['passed'] and o['status']=='terminal' and o['child_exit']==0 and o['source_unchanged'] and o['whole_inputs_unchanged'] and not o['remaining_owned_processes'] and not o['cleanup_errors'] and all(x['absent'] for x in o['observed_births'])
  mapping.extend([('out/forensic-alcubierre-png','files/out/forensic-alcubierre-png'),('out/forensic-alcubierre-completion','files/out/forensic-alcubierre-completion')])
 message='Preserve full-size reconstructed legacy and ported '+label+' identity (#100)\n'
else:raise SystemExit(kind)
assert not bundle.exists();bundle.mkdir(parents=True)
for a,b in mapping:
 p=root/a;q=bundle/b;q.parent.mkdir(parents=True,exist_ok=True)
 if p.is_dir():shutil.copytree(p,q,symlinks=True)
 else:shutil.copy2(p,q)
shutil.copy2(Path(__file__),bundle/'preserve-remediation-bundle.py')
for extra in ('out/capture-forensic-live.py','out/verify-preserved-remediation.py'):
 p=root/extra
 if p.is_file():
  q=bundle/'helpers'/p.name;q.parent.mkdir(exist_ok=True);shutil.copy2(p,q)
write(bundle/'result.json',summary)
files=[p for p in sorted(bundle.rglob('*')) if p.is_file() or p.is_symlink()]
manifest={'schema':'sirius-complete-raw-artifact-bundle-v1','qualification_claimed':False,'files':{str(p.relative_to(bundle)):identity(p) for p in files}}
write(bundle/'manifest.json',manifest)
fd,index=tempfile.mkstemp(prefix='preserve-'+kind+'-',dir=root/'out');os.close(fd);os.unlink(index);env=dict(os.environ,GIT_INDEX_FILE=index)
try:
 subprocess.run(['git','read-tree','--empty'],env=env,check=True);rows=[]
 for p in sorted(bundle.rglob('*')):
  if p.is_symlink():mode='120000';sha=subprocess.check_output(['git','hash-object','-w','--stdin'],input=os.readlink(p).encode()).decode().strip()
  elif p.is_file():mode='100755' if p.stat().st_mode&0o111 else '100644';sha=subprocess.check_output(['git','hash-object','-w','--',str(p)]).decode().strip()
  else:continue
  rows.append(f'{mode} {sha}\t{p.relative_to(root)}\0')
 subprocess.run(['git','update-index','-z','--index-info'],input=''.join(rows).encode(),env=env,check=True);tree=subprocess.check_output(['git','write-tree'],env=env).decode().strip();commit=subprocess.check_output(['git','commit-tree',tree,'-p',source],input=message.encode(),env=env).decode().strip();subprocess.run(['git','tag',tag,commit],check=True)
 result={'bundle':str(bundle.relative_to(root)),'tag':tag,'commit':commit,'files':len(rows),'manifest':identity(bundle/'manifest.json')}
 if kind=='translation':result['provider_payload_tree']=subprocess.check_output(['git','rev-parse',commit+':'+str(bundle.relative_to(root))+'/files/out/macos-moltenvk-unforced']).decode().strip()
 for p in sorted(bundle.rglob('*')):
  if p.is_symlink():expected=os.readlink(p).encode()
  elif p.is_file():expected=p.read_bytes()
  else:continue
  assert subprocess.check_output(['git','show',commit+':'+str(p.relative_to(root))])==expected
 result['all_git_bytes_verified']=True
 print(json.dumps(result,indent=2))
finally:
 if Path(index).exists():Path(index).unlink()
