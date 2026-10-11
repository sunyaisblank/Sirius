"""Freeze the original twenty-condition default ABBA gate before its first invocation."""
from pathlib import Path
import hashlib,json,subprocess,sys
sys.dont_write_bytecode=True
W=Path(__file__).resolve().parent;R=W.parents[1]
HEAD='b81e061fa081cdaed37eca3add950bd254bf296d'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def dump(p,d):Path(p).write_text(json.dumps(d,indent=2)+'\n')
def identity(p):return {'bytes':Path(p).stat().st_size,'sha256':sha(p)}
assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()==HEAD
assert not subprocess.check_output(['git','status','--porcelain'],cwd=R)
assert not (W/'retention-gate-frozen.json').exists()
assert not any((W/n).exists() for n in ['a1','b1','b2','a2'])
controls=['candidate-cache-primer','candidate-mixed-cached','candidate-rk-cached','candidate-schwarzschild-cached']
additional=[W/n for n in ['prepare_gate.py','compare_software.py','run_sequence.py','read_controls.py','run_cached_owned.py','numerical-control-readback.json','source-independent-review.json','compiled-independent-review.json','host-apparatus-independent-review.json','cold-rk-independent-review.json']]
for name in controls:
    folder=W/name;owner=load(folder/'owner.json');reader=load(folder/'control-readback.json')
    assert owner['passed'] and owner['child_exit']==0 and owner['stop_reason'] is None
    assert owner['remaining_owned_processes']==[] and owner['cleanup_errors']==[] and all(r['absent'] for r in owner['observed_births'])
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert reader['pass'] and reader['owner']==identity(folder/'owner.json') and reader['xml']==identity(folder/'gtest.xml')
    before,after=load(folder/'inputs-before.json'),load(folder/'inputs-after.json');assert before==after
    for key,row in before.items():
        path=W/'action-versions'/(name+'.json') if key==str((W/'actions.json').relative_to(R)) else Path(row['resolved_path'])
        assert identity(path)=={k:row[k] for k in ('bytes','sha256')},str(path)
    additional += [folder/n for n in ['owner.json','control-readback.json','gtest.xml','inputs-before.json','inputs-after.json']]
host=load(W/'host-preparation.json');assert host['pass'] and host['changed_arrays']==['kTransportPortableNormalSumShader'] and host['unchanged_other_arrays']==39
original=R/'attestations/software-vulkan/dc714b9/guarded-upmultiply-trial/diagnostic/retention-gate-frozen.json'
prior=load(original)
assert prior['sequence']==['a1','b1','b2','a2'] and prior['within_slot_repeats']==3
assert prior['complete_call_minimum_improvement_percent']==5 and prior['every_inclusive_stage_minimum_improvement_percent']==-5 and prior['all20_gates_required']
additional += [original,W/'retention-gate-frozen.json']
actions=load(W/'actions.json');template=actions['candidate-timestamps']
executables={}
for producer in ['baseline','candidate']:
    path=W/(producer+'-backend-tests');record=host['artifacts'][producer]['executable']
    assert identity(path)=={k:record[k] for k in ('bytes','sha256')}
    executables[producer]={'path':str(path.relative_to(R)),'source_revision':HEAD,**identity(path)}
for label in ['a1','b1','b2','a2']:
    producer='baseline' if label.startswith('a') else 'candidate'
    action=json.loads(json.dumps(template));action['argv']=[executables[producer]['path'],'--gtest_filter=RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults']
    action['scope']='One frozen '+label+' slot of the original default ABBA twenty-condition software benefit gate, cache disabled, exact original observer/scientific assertions/32 complete buffers/51 ordered markers. No adoption/rendering/final qualification.'
    action['frozen_inputs']=sorted(set(action['frozen_inputs']+[str(p.resolve()) for p in additional]))
    actions[label]=action
dump(W/'actions.json',actions)
for label in ['a1','b1','b2','a2']:(W/'action-versions'/(label+'.json')).write_bytes((W/'actions.json').read_bytes())
freeze={'candidate_revision':HEAD,'baseline_revision':HEAD,'original_baseline_producer_revision':'7f028007a4935ffdefe5ee9c886b138299f50dd8',
    'candidate_kind':'Unadopted compact full-table Transport evaluator in a private factory host; production source/head unchanged',
    'candidate_source':identity(W/'kernels/retained_transport.slang'),'candidate_module':identity(W/'literal-transport.spv'),
    'executables':executables,'sequence':['a1','b1','b2','a2'],'within_slot_repeats':3,
    'complete_call_minimum_improvement_percent':5,'every_inclusive_stage_minimum_improvement_percent':-5,'all20_gates_required':True,
    'reader_sha256':sha(W/'software_contracts.py'),'comparator_sha256':sha(W/'compare_software.py'),'launcher_sha256':sha(W/'run_sequence.py'),
    'actions_sha256':sha(W/'actions.json'),'owner_sha256':sha(W/'run_owned.py'),'guard_sha256':sha(W/'owned_guard.py'),'preparer_sha256':sha(__file__),
    'prior_original_gate':{'path':str(original.relative_to(R)),**identity(original)},'changed_selected_modules':['transport'],
    'scope':'Frozen before any benefit invocation. Original four complete-call >=5% and sixteen inclusive-stage >=-5% conditions; fixed forward/reverse A1/B1/B2/A2 and three original repeats. Actual baseline/candidate same strict private factory frames, only Transport differs; cold cache-disabled ownership remains separate from prepared-cache correctness. No retry pooling, scene/threshold/assertion changes, product adoption, rendering or final qualification.'}
dump(W/'retention-gate-frozen.json',freeze)
print(json.dumps({'gate_sha256':sha(W/'retention-gate-frozen.json'),'actions_sha256':sha(W/'actions.json'),'sequence':freeze['sequence'],'conditions':20,'product_adopted':False}))
