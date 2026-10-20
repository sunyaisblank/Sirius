"""Preserve exact issue readbacks and verified scratch retirement on the evidence branch."""
from pathlib import Path
import hashlib,json,os,shutil,stat,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;SOURCE='630844e775a85e74eb5d9ed328a2b8019abb640d';PARENT='121db07cabe4508890c3ba3045bfe040805df8be';BRANCH='development/software-single-invocation-2026-10-20';T=R/'attestations/qualification/630844e7/software-single-invocation-cost-closeout';INDEX=W/'publication.index'
def git(*a,data=None,env=None):return subprocess.check_output(['git',*a],input=data,env=env)
def identity(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest(),'mode':oct(stat.S_IMODE(p.stat().st_mode))}
assert git('rev-parse','HEAD').decode().strip()==SOURCE and git('status','--porcelain')==b''
assert git('rev-parse',BRANCH).decode().strip()==PARENT and not T.exists()
v=json.loads((W/'publication-verification.json').read_text());r=json.loads((W/'retirement.json').read_text())
assert v['comment_exact_readback'] and v['issue44_body_and_original_acceptance_unchanged'] and v['PR62_open_draft_source_unchanged']
assert r['passed'] and r['retired_files']==402 and all(not (R/n).exists() for n in r['retired_roots'])
state={'source_branch':git('branch','--show-current').decode().strip(),'source_head':SOURCE,'source_status':git('status','--porcelain').decode(),'source_tracking_counts':git('rev-list','--left-right','--count','HEAD...@{upstream}').decode().strip(),'worktrees':git('worktree','list','--porcelain').decode(),'prototype_head_before_closeout':PARENT,'all_original_scientific_runtime_and_release_gates_remain_required':True,'goal_status':'active'}
assert state['source_tracking_counts']=='0\t0' and state['worktrees'].count('worktree ')==1
(W/'source-workspace-before-closeout.json').write_text(json.dumps(state,indent=2)+'\n')
files={p.name:identity(p) for p in W.iterdir() if p.is_file() and p!=INDEX};T.mkdir(parents=True)
for n,pin in files.items():shutil.copy2(W/n,T/n);assert identity(T/n)==pin
(T/'MANIFEST.json').write_text(json.dumps({'schema':'sirius-single-invocation-cost-publication-closeout-v1','operator_source_revision':SOURCE,'parent_evidence_checkpoint':PARENT,'files':files,'source_root':str(W.relative_to(R)),'publication':v,'retirement':r,'prototype_rejected':True,'source_changes':False,'production_adoption':False,'full_qualification':False,'current_goal_status':'active'},indent=2)+'\n')
env=dict(os.environ,GIT_INDEX_FILE=str(INDEX))
try:
 git('read-tree',PARENT,env=env);rows=[]
 for p in sorted(T.iterdir()):
  blob=git('hash-object','-w',str(p)).decode().strip();rows.append(('100644 '+blob+'\t'+str(p.relative_to(R))+'\n').encode())
 git('update-index','--index-info',data=b''.join(rows),env=env);tree=git('write-tree',env=env).decode().strip();c=git('commit-tree',tree,'-p',PARENT,data=b'Close Endpoint cost screen with exact issue readback and Git-backed retirement\n').decode().strip();git('update-ref','refs/heads/'+BRANCH,c,PARENT)
 for p in T.iterdir():
  b=git('show',c+':'+str(p.relative_to(R)));assert len(b)==p.stat().st_size and hashlib.sha256(b).hexdigest()==identity(p)['sha256']
finally:INDEX.unlink(missing_ok=True)
assert all(identity(W/n)==pin for n,pin in files.items()) and git('status','--porcelain')==b''
print(json.dumps({'commit':c,'branch':BRANCH,'bundle':str(T.relative_to(R)),'files':len(files),'issue_comment':v['comment_url'],'retired_prototype_files':r['retired_files'],'source_clean_and_pushed':True}))
