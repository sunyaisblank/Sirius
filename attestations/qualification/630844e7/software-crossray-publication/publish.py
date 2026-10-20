"""Preserve exact queue publication and verified predecessor retirement."""
from pathlib import Path
import hashlib,json,os,shutil,stat,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;SOURCE='630844e775a85e74eb5d9ed328a2b8019abb640d';PARENT='65a466249d923a4b8c2ce1a02496cafdaa2e733a';BR='development/software-crossray-endpoint-2026-10-20';T=R/'attestations/qualification/630844e7/software-crossray-publication';INDEX=W/'publication.index'
def git(*a,data=None,env=None):return subprocess.check_output(['git',*a],input=data,env=env)
def signature(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest(),'mode':oct(stat.S_IMODE(p.stat().st_mode))}
assert git('rev-parse','HEAD').decode().strip()==SOURCE and git('status','--porcelain')==b'' and git('rev-parse',BR).decode().strip()==PARENT and not T.exists()
v=json.loads((W/'publication-verification.json').read_text());r=json.loads((W/'retirement.json').read_text());assert v['comment_exact_readback'] and v['issue44_original_body_and_acceptance_unchanged'] and v['PR62_open_draft_source_unchanged'] and r['passed']
assert r['retired_files']==192 and all(not (R/p).exists() for p in r['retired_roots']) and (R/r['corrected_root']).is_dir()
state={'source_branch':git('branch','--show-current').decode().strip(),'source_head':SOURCE,'source_status':git('status','--porcelain').decode(),'source_tracking_counts':git('rev-list','--left-right','--count','HEAD...@{upstream}').decode().strip(),'prototype_branch':BR,'prototype_head_before_publication':PARENT,'prototype_tracking_counts_before_publication':git('rev-list','--left-right','--count',BR+'...origin/'+BR).decode().strip(),'worktrees':git('worktree','list','--porcelain').decode(),'retained_corrected_root':r['corrected_root'],'retained_corrected_files':r['corrected_root_files_retained'],'goal_status':'active','all_original_scientific_runtime_frame_and_installed_release_obligations_remain_required':True}
assert state['source_tracking_counts']==state['prototype_tracking_counts_before_publication']=='0\t0' and state['worktrees'].count('worktree ')==1
(W/'source-workspace-before-publication.json').write_text(json.dumps(state,indent=2)+'\n')
files={p.name:signature(p) for p in W.iterdir() if p.is_file() and p!=INDEX};T.mkdir(parents=True)
for n,pin in files.items():shutil.copy2(W/n,T/n);assert signature(T/n)==pin
(T/'MANIFEST.json').write_text(json.dumps({'schema':'sirius-crossray-queue-publication-retirement-v1','operator_source_revision':SOURCE,'parent_corrected_source_checkpoint':PARENT,'files':files,'source_root':str(W.relative_to(R)),'publication':v,'retirement':r,'first_variant_disposition':'REJECT_EIGHT_PRIVATE_RAY_ENDPOINT_ACTIVE_BOUND_AND_COST','corrected_variant_compiler_and_emitted_guard_only':True,'corrected_tail_observer_source_prepared':True,'corrected_tail_observer_compiled_or_run':False,'production_adopted':False,'full_qualification':False,'goal_status':'active'},indent=2)+'\n')
env=dict(os.environ,GIT_INDEX_FILE=str(INDEX))
try:
 git('read-tree',PARENT,env=env);rows=[]
 for p in sorted(T.iterdir()):
  blob=git('hash-object','-w',str(p)).decode().strip();rows.append(('100644 '+blob+'\t'+str(p.relative_to(R))+'\n').encode())
 git('update-index','--index-info',data=b''.join(rows),env=env);tree=git('write-tree',env=env).decode().strip();c=git('commit-tree',tree,'-p',PARENT,data=b'Publish cross-ray ownership correction and retire rejected scoped workspace\n').decode().strip();git('update-ref','refs/heads/'+BR,c,PARENT)
 for p in T.iterdir():
  data=git('show',c+':'+str(p.relative_to(R)));pin=signature(p);assert len(data)==pin['bytes'] and hashlib.sha256(data).hexdigest()==pin['sha256']
finally:INDEX.unlink(missing_ok=True)
assert all(signature(W/n)==pin for n,pin in files.items()) and git('status','--porcelain')==b''
print(json.dumps({'commit':c,'branch':BR,'bundle':str(T.relative_to(R)),'files':len(files),'issue_comment':v['comment_url'],'rejected_files_retired':r['retired_files'],'corrected_source_files_retained':r['corrected_root_files_retained'],'source_clean_and_pushed':True}))
