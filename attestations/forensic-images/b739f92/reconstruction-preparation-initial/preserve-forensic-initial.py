import hashlib,json,os,pathlib,shutil,subprocess,tempfile
root=pathlib.Path.cwd();work=root/'out/forensic-reconstruction-b739f92';bundle=root/'attestations/forensic-images/b739f92/reconstruction-preparation-initial'
assert not bundle.exists()
def identity(p):
 b=p.read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
assert identity(work/'preparation.json')['sha256']=='931a6cb3b513e6afd355c214397e389b91ad60f52983a9f5ae16c8cfca1a7e35'
assert identity(work/'pins.json')['sha256']=='0e4db79d3787b780a64e19be6b1be8c5aa5457f5c10aa7663e236ab4a3a33214'
for name,code,passed in [('build-comparator',0,True),('comparator-controls',0,True),('configure-renderers',1,False)]:
 o=json.loads((work/name/'owner.json').read_text());assert o['status']=='terminal' and o['child_exit']==code and o['passed']==passed and not o['remaining_owned_processes'] and not o['cleanup_errors'] and o['whole_inputs_unchanged'] and o['source_unchanged']
 assert o['observed_births'] and all(x['absent'] for x in o['observed_births'])
summary={'schema':'sirius-forensic-reconstruction-preparation-v1','evidence_kind':'NEW_RECONSTRUCTION','runner_revision':'5733117ba8ba4430d723d62294f991db423fbf4c','source_export_revision':'b739f920c2ae67847ed7961e011d2350a839e410','legacy_runtime_origin':'ec2d8e506dfd0477e3c9eb4b2b4aa1bee9831fe2','historical_tape_identity_claimed':False,'renders_executed':False,'comparator_built':True,'comparator_comparisons_passed':17,'comparator_fixture_creation_passed':True,'comparator_domains':['PNG fullRGBA8','EXR decoded IEEE754 binary32 words from four ZIP HALF channels','PPM fullP6 RGB8'],'renderer_configuration':'FAILED missing historical lib/glfw/CMake/glfw3Config.cmake.in; distinct build-scaffold repair pending','renderers_built':False,'required_full_scene_rays_not_executed':291635200,'qualification_claimed':False}
bundle.mkdir(parents=True);shutil.copytree(work,bundle/'packet');shutil.copytree(root/'bin/linux-gcc/forensic-reconstruction-b739f92',bundle/'build-products')
shutil.copy2(root/'out/preserve-forensic-initial.py',bundle/'preserve-forensic-initial.py')
(bundle/'result.json').write_text(json.dumps(summary,indent=2)+'\n')
manifest={'schema':'sirius-raw-reconstruction-bundle-v1','qualification_claimed':False,'files':{str(p.relative_to(bundle)):identity(p) for p in sorted(bundle.rglob('*')) if p.is_file()}}
(bundle/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
fd,index=tempfile.mkstemp(prefix='forensic-initial-',dir=root/'out');os.close(fd);os.unlink(index);env=dict(os.environ,GIT_INDEX_FILE=index)
try:
 subprocess.run(['git','read-tree','--empty'],env=env,check=True);rows=[]
 for p in sorted(bundle.rglob('*')):
  if not p.is_file():continue
  assert not p.is_symlink(),p
  sha=subprocess.check_output(['git','hash-object','-w','--',str(p)]).decode().strip();mode='100755' if p.stat().st_mode&0o111 else '100644';rows.append(f'{mode} {sha}\t{p.relative_to(root)}\0')
 subprocess.run(['git','update-index','-z','--index-info'],input=''.join(rows).encode(),env=env,check=True)
 tree=subprocess.check_output(['git','write-tree'],env=env).decode().strip();commit=subprocess.check_output(['git','commit-tree',tree,'-p',summary['runner_revision']],input=b'Preserve historical reconstruction comparator controls and missing scaffold failure (#100)\n',env=env).decode().strip()
 tag='evidence/issue-100-reconstruction-preparation-initial';subprocess.run(['git','tag',tag,commit],check=True)
 print(json.dumps({'tag':tag,'commit':commit,'files':len(rows),'manifest':identity(bundle/'manifest.json'),'result':summary},indent=2))
finally:
 if pathlib.Path(index).exists():pathlib.Path(index).unlink()
