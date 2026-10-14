from pathlib import Path
import hashlib,json,os,subprocess,tempfile
root=Path.cwd();w=root/'out/native-provider-integration';m=json.loads((w/'candidate.json').read_text());checks=json.loads((w/'root-integrated-source-checks.json').read_text());assert checks['candidate_inputs']==m['files'];assert checks['actual_candidate_workflow_policy'] and checks['new12_source_boundary_controls'];parent=m['base_revision'];branch='refs/heads/development/native-macos-pinned-provider-2026-10-14'
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==m['checkout_unchanged_revision'];assert not subprocess.check_output(['git','status','--porcelain']).strip()
review=json.loads((w/'independent-workflow-review.json').read_text());assert review['material_blockers']==[]
fd,index=tempfile.mkstemp(dir=root/'out',prefix='native-provider-index-');os.close(fd);os.unlink(index);env=dict(os.environ,GIT_INDEX_FILE=index)
try:
 subprocess.run(['git','read-tree',parent],env=env,check=True)
 for name,row in m['files'].items():
  p=w/name;b=p.read_bytes();assert {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}==row['candidate'];blob=subprocess.check_output(['git','hash-object','-w','--',str(p)],text=True).strip();subprocess.run(['git','update-index','--add','--cacheinfo',row['mode'],blob,name],env=env,check=True)
 tree=subprocess.check_output(['git','write-tree'],env=env,text=True).strip();assert set(subprocess.check_output(['git','diff','--name-only',parent,tree],text=True).splitlines())==set(m['files']);subprocess.run(['git','diff','--check',parent,tree],check=True)
 product_diff=subprocess.check_output(['git','diff','--name-only',m['product_source'],tree],text=True).splitlines();assert set(product_diff)==set(m['diagnostic_source_whitelist'])
 commit=subprocess.check_output(['git','commit-tree',tree,'-p',parent],env=env,input=b'Provision and bind pinned MoltenVK in current-source native qualification (#48)\n').decode().strip()
 for name in m['files']:assert subprocess.check_output(['git','show',commit+':'+name])==(w/name).read_bytes()
 subprocess.run(['git','update-ref',branch,commit,'0'*40],check=True)
 result={'commit':commit,'parent':parent,'branch':branch,'working_checkout_unchanged':True,'checks':checks,'independent_component_review':json.loads((w/'independent-component-review.json').read_text()),'independent_workflow_review':review,'product_source_difference_exactly9reviewed_nonproduct_paths':product_diff,'nativeexecution_or_qualification_claimed':False};(w/'source-commit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'commit':commit,'branch':branch,'changed_files':list(m['files']),'product_difference':product_diff}))
finally:
 if Path(index).exists():Path(index).unlink()
