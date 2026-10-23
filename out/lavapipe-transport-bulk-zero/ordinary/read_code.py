"""Reuse the complete strict code gate for the exact Transport scratch array."""
from pathlib import Path
import hashlib
_basis=Path.cwd()/'out/lavapipe-owner-bulk-zero/ordinary/read_code.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest()=='5eb4d91076d7566965ce3734d593eb8c224e8f06cda3306ace3dfcaa371a3978'
_text=_basis.read_text()
_prefix=_text.split('\ncompiler=',1)[0].replace('%reg','%reg32')
for _old,_new in {
 '%.fca.{i}.gep = getelementptr':'%.fca.{i}.gep{12990+i if i < 125 else ""} = getelementptr',
 'ptr %.fca.{i}.gep, align 32':'ptr %.fca.{i}.gep{12990+i if i < 125 else ""}, align 32',
}.items():
 assert _prefix.count(_old)==1;_prefix=_prefix.replace(_old,_new)
_start="for label,kind,count in [('endpoint','Endpoint',3060),('dopri_phase','DopriPhase',3145)]:"
_loop=_text[_text.index(_start):_text.index('\nout=')]
_loop=_loop.replace(_start,"for label,kind,count in [('transport','Transport',2295)]:")
_old="    assert current_norm==prior_norm, ('IR changed outside declared full-clear equivalence',kind)"
_new="""    if current_norm!=prior_norm:
        import difflib
        (W/'transport-remaining-IR.diff').write_text(''.join(difflib.unified_diff(prior_norm.splitlines(True),current_norm.splitlines(True),fromfile='prior',tofile='candidate')))
        (W/'strict-equivalence-refusal.json').write_text(json.dumps({'passed':False,'actual_complete_clear_observed':True,'actual_clear_bytes':clear['bytes'],'reason':'Remaining optimized IR differs; passive source/IR review required before any numerical run.','no_capture_or_workload_retry_authorized':True},indent=2)+'\\n')
    assert current_norm==prior_norm, ('IR changed outside declared full-clear equivalence',kind)"""
assert _loop.count(_old)==1;_loop=_loop.replace(_old,_new)
_setup="""
compiler=doc(R/'out/native-cpu-row-ownership/compiler-result.json')
provider=doc(W.parent/'provider-build-result.json');assert provider['passed']
driver=R/provider['library'];assert pin(driver)==provider['library_identity']
rows={}
"""
_end="""
out={'structural_premise_observed':True,'rows':rows,'complete_bulk_clear_observed':True,
     'remaining_optimized_IR_equivalent':True,'original_Default_machine_pipeline_observed':True,
     'actual_provider':{'path':str(driver.relative_to(R)),**pin(driver)},
     'scope':'One full row-T4B guard/canary. Complete clear observed; no retained numerical/benefit/frame/science/release claim. E/DP source/code proof remains separate.'}
target=W/'ordinary-code-review.json';assert not target.exists()
target.write_text(json.dumps(out,indent=2)+'\\n')
print(json.dumps({'passed':True,'Transport_clear_bytes':73440}))
"""
exec(compile(_prefix+'\n'+_setup+'\n'+_loop+'\n'+_end,str(_basis),'exec'))
