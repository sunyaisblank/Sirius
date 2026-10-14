"""Seal one of three fixed completed cleanup receipts using an isolated Git index."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,tempfile
ROOT=Path.cwd()
SOURCE='4498cc41168ad29fc4c111819a4cb7ccaa54a8bf'
KINDS={
 'fp64-estimate-terminal-scratch':[
  'out/portable-divide-fp64-estimate-screen','bin/linux-gcc/portable-divide-fp64-estimate-screen',
  'out/portable-divide-fp64-estimate-property-query'],
 'native-terminal-scratch':[
  'out/native-windows-system-5b5fad9','bin/windows-msvc/native-system-5b5fad9',
  'out/macos-current-provider-recovery','out/native-provider-integration/full-ci-38577826600/linux-gate',
  'out/native-provider-integration/full-ci-38577826600/macos-build'],
 'minimum-step-floor-terminal-scratch':[
  'out/minimum-step-floor-current','out/minimum-step-floor-registration','out/minimum-step-ellis-floor',
  'out/macos-failure-diagnostics','out/minimum-step-floor-repair','out/minimum-step-floor-verification']}
kind=sys.argv[1]
assert kind in KINDS
BUNDLE=ROOT/'attestations/cleanup/2026-10-15'/kind
TAG='evidence/issue-'+('44-' if kind=='fp64-estimate-terminal-scratch' else '48-' if kind=='native-terminal-scratch' else '53-')+kind+'-cleanup'
def git(*args,env=None):return subprocess.check_output(['git','-c','gc.auto=0',*args],env=env)
def pin(path):
 raw=path.read_bytes();return {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
assert git('rev-parse','HEAD').decode().strip()==SOURCE and not git('status','--porcelain').strip()
r=json.loads((BUNDLE/'receipt.json').read_text())
assert r['status']=='preserved-terminal-roots-removed' and r['source_unchanged'] is True
recorded_source=r.get('source',r.get('runner_revision'))
assert (recorded_source.get('revision') if isinstance(recorded_source,dict) else recorded_source)==SOURCE
assert r['qualification_claimed'] is False and r['global_consumer_absence_claimed'] is False
removed=[x if isinstance(x,str) else x['path'] for x in r['removed']]
assert removed==KINDS[kind]
assert all(not (ROOT/name).exists() for name in removed)
assert (ROOT/'bin/linux-gcc').is_dir() and (ROOT/'bin/linux-gcc-sanitize').is_dir()
assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+TAG]).returncode==1
if kind=='fp64-estimate-terminal-scratch':
 apparatus=ROOT/'out/portable-divide-fp64-estimate-retirement'
 assert r['recorded_owner_count']==5 and r['recorded_birth_count']==9
 assert git('rev-list','--left-right','--count','HEAD...@{upstream}').split()==[b'0',b'0']
 assert [v for v in git('worktree','list','--porcelain').decode().splitlines() if v.startswith('worktree ')]==['worktree '+str(ROOT)]
 for name in ('independent-source-review.json','independent-audit-review.json'):
  assert json.loads((apparatus/name).read_text())['passed'] is True
 assert not (BUNDLE/'apparatus').exists()
 shutil.copytree(apparatus,BUNDLE/'apparatus',symlinks=True)
shutil.copy2(Path(__file__),BUNDLE/Path(__file__).name)
files=[p for p in sorted(BUNDLE.rglob('*')) if p.is_file() and p!=BUNDLE/'manifest.json']
assert files and all(not p.is_symlink() for p in files)
(BUNDLE/'manifest.json').write_text(json.dumps({'qualification_claimed':False,'files':{str(p.relative_to(BUNDLE)):pin(p) for p in files}},indent=2,sort_keys=True)+'\n')
files.append(BUNDLE/'manifest.json')
fd,index=tempfile.mkstemp(dir=ROOT/'out',prefix='fixed-cleanup-index-');os.close(fd);os.unlink(index)
env=dict(os.environ,GIT_INDEX_FILE=index)
try:
 git('read-tree','--empty',env=env)
 rows=[]
 for p in sorted(files):
  blob=git('hash-object','-w','--',str(p)).decode().strip()
  rows.append(f'100644 {blob}\t{p.relative_to(ROOT)}\0')
 subprocess.run(['git','update-index','-z','--index-info'],input=''.join(rows).encode(),env=env,check=True)
 tree=git('write-tree',env=env).decode().strip()
 commit=subprocess.check_output(['git','-c','gc.auto=0','commit-tree',tree,'-p',SOURCE],input=('Preserve completed '+kind+' cleanup\n').encode(),env=env).decode().strip()
 subprocess.run(['git','tag',TAG,commit],check=True)
 for p in files:assert git('show',commit+':'+str(p.relative_to(ROOT)))==p.read_bytes()
 assert not git('status','--porcelain').strip()
 print(json.dumps({'tag':TAG,'commit':commit,'bundle':str(BUNDLE.relative_to(ROOT)),'files':len(files),'manifest':pin(BUNDLE/'manifest.json'),'all_git_bytes_verified':True,'pushed':False}))
finally:
 if Path(index).exists():Path(index).unlink()
