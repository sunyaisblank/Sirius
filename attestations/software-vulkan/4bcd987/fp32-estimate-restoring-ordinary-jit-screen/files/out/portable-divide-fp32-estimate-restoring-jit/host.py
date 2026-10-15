"""Reuse canonical host/stock; admit one full FP32-restoring-estimate ordinary canary capture."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys
sys.dont_write_bytecode=True
from prerequisites import WORK,ROOT,P,load,record,content,module,accepted,verify_graph,verify_stock
BUILD=ROOT/'bin/linux-gcc/portable-divide-fp32-estimate-restoring-jit'

def write(path,value):
    path=Path(path);assert not path.exists(),path
    path.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')

def canonical_json(value):
    # Exact durable JSON normalization: integer object keys only; every field/type/array order retained.
    normalized=json.loads(json.dumps(value,allow_nan=False))
    return json.dumps(normalized,sort_keys=True,separators=(',',':'),allow_nan=False)

def verify_static():
    for name,pin in load(WORK/'source-freeze.json')['pins'].items():assert record(name)==pin,name

def verify_providers():
    for name,pin in P['providers'].items():
        path=Path(name);current=path.stat()
        value=record(path);assert {k:value[k] for k in ('bytes','sha256','resolved_path')}=={k:pin[k] for k in ('bytes','sha256','resolved_path')}
        assert [os.major(current.st_dev),os.minor(current.st_dev)]==pin['device'] and current.st_ino==pin['inode']

def fixed_context():
    context=load(WORK/'runtime-context.json');assert context['revision']==P['operator_revision']
    assert context['source_freeze']==record(WORK/'source-freeze.json') and context['preparation']==record(WORK/'preparation.json')
    assert context['capture_actions']==record(WORK/'capture-actions.json') and context['reuse_acceptance']==record(WORK/'reuse-acceptance.json')
    assert context['host']==record(BUILD/'short_input_jit_probe') and context['modules']['candidate']==record(WORK/'module-inputs/candidate-endpoint-normal.spv')
    actions=load(WORK/'capture-actions.json');expected=load(WORK/'capture-actions.unpinned.json')
    assert expected['capture-candidate-endpoint-normal']['revision']=='PENDING_OPERATOR_CLEAN_CONTEXT'
    expected['capture-candidate-endpoint-normal']['revision']=P['operator_revision'];assert actions==expected
    return context

def pin_runtime():
    verify_static();stock=verify_stock();graph=verify_graph();verify_providers()
    assert not BUILD.exists() and not (WORK/'module-inputs').exists()
    assert not any((WORK/name).exists() for name in ('runtime-context.json','capture-actions.json','reuse-acceptance.json','capture-candidate-endpoint-normal','capture-acceptance.json'))
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()==P['operator_revision']
    assert not subprocess.check_output(['git','status','--porcelain=v1'],cwd=ROOT)
    BUILD.mkdir();executable=BUILD/'short_input_jit_probe';shutil.copyfile(P['canary'],executable);executable.chmod(0o755)
    assert content(executable)==P['canary_identity']
    (WORK/'module-inputs').mkdir();candidate=WORK/'module-inputs/candidate-endpoint-normal.spv'
    shutil.copyfile(graph['candidate_module']['resolved_path'],candidate);assert content(candidate)==P['candidate_module']
    actions=load(WORK/'capture-actions.unpinned.json');assert list(actions)==['capture-candidate-endpoint-normal']
    assert actions['capture-candidate-endpoint-normal']['revision']=='PENDING_OPERATOR_CLEAN_CONTEXT'
    actions['capture-candidate-endpoint-normal']['revision']=P['operator_revision'];write(WORK/'capture-actions.json',actions)
    write(WORK/'reuse-acceptance.json',{'passed':True,'stock':stock,'graph':graph,'host_recompiled':False,'stock_recaptured':False,
          'feature_scope':'Unmodified canonical host: Int64 disabled; stock default Float64 enable remains unused by this only32-bit shader. New actual arguments must be observed; old canonical stock has no raw55-argument capture.'})
    write(WORK/'runtime-context.json',{'revision':P['operator_revision'],'host':record(executable),'modules':{'candidate':record(candidate)},
          'capture_actions':record(WORK/'capture-actions.json'),'source_freeze':record(WORK/'source-freeze.json'),'preparation':record(WORK/'preparation.json'),
          'reuse_acceptance':record(WORK/'reuse-acceptance.json'),
          'scope':'One original early-return short-input candidate capture. No retained arithmetic, Float32 accuracy/hit rate, whole-call performance or qualification claim.'})
    print('Exact canonical host/module runtime inputs pinned; no provider process launched.',flush=True)

def verify_new_maps(folder,owner,code):
    verify_providers();pid=load(folder/'inferior.json')['pid'];births={row['pid']:row['saved_start_ticks'] for row in owner['observed_births']}
    assert pid in births;expected=set(P['providers'])
    def check(lines):
        found=set()
        for line in lines:
            fields=line.split(None,5)
            if len(fields)!=6 or fields[5] not in expected:continue
            pin=P['providers'][fields[5]]
            assert [int(value,16) for value in fields[3].split(':')]==pin['device'] and int(fields[4])==pin['inode'];found.add(fields[5])
        return found
    observations=[row for row in owner['provider_mappings'] if row['pid']==pid]
    assert observations and all(row['start_ticks']==births[pid] and row['owned_session']==owner['child_pid'] for row in observations)
    assert set().union(*(check(row['maps']) for row in observations))==expected
    events=[row for row in load(folder/'bindings.json')['jit_events'] if row.get('file')==code['actual_entry']['file']]
    assert len(events)==1 and check(events[0]['maps_at_first_call'].splitlines())==expected

def review():
    verify_static();stock=verify_stock();graph=verify_graph();verify_providers();context=fixed_context()
    reuse=load(WORK/'reuse-acceptance.json');assert reuse['stock']==stock
    # Fresh /proc observations are rechecked, not compared as a stable historical identity.
    assert {k:v for k,v in reuse['graph'].items() if k!='fresh_graph_birth_absence'}=={k:v for k,v in graph.items() if k!='fresh_graph_birth_absence'}
    assert all(row['saved_birth_absent'] for row in reuse['graph']['fresh_graph_birth_absence'])
    folder=WORK/'capture-candidate-endpoint-normal';owner=accepted(folder)
    action=load(WORK/'capture-actions.json')['capture-candidate-endpoint-normal']
    assert owner['source']=={'revision':P['operator_revision'],'status':''} and owner['limits']==P['limits']==action['limits']
    assert owner['argv']==[arg.replace('{folder}',str(folder)) for arg in action['argv']]
    before=load(folder/'inputs-before.json')
    for name,pin in before.items():
        path=Path(name);path=path if path.is_absolute() else ROOT/path
        value=record(path);assert {k:value[k] for k in pin}==pin,name
    for name in action['frozen_inputs']:
        path=Path(name);key=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
        value=record(path);assert {k:value[k] for k in before[key]}==before[key]
    analyzer=module('unchanged_ordinary_code_reader',WORK/'analyze_capture.py')
    feature_reader=module('canonical_actual_feature_reader',WORK/'read_features.py')
    code=analyzer.inspect(folder);features=feature_reader.inspect(folder)
    assert canonical_json(features)==canonical_json(load(folder/'feature-verification.json'))
    stock_folder=Path(stock['stock_capture_folder']);stock_code=load(stock_folder/'actual-code-review.json')
    assert canonical_json(analyzer.inspect(stock_folder))==canonical_json(stock_code)
    assert code['ordinary_pipelines']==P['ordinary_pipelines']==stock_code['ordinary_pipelines']
    assert len(features['features'])==len(load(WORK/'feature-layout.json')['features_fields'])==55
    assert features['features']['shaderInt64']==0 and features['features']['shaderFloat64']==1
    assert features['actual_public_create_arguments_enable_shaderInt64'] is False
    assert features['successful_logical_device_used_for_actual_complete_shader']
    assert code['actual_loaded_module']['sha256']==P['candidate_module']['sha256'] and code['actual_loaded_module']['bytes']==P['candidate_module']['bytes']
    assert content(BUILD/'short_input_jit_probe')==P['canary_identity']
    verify_new_maps(folder,owner,code)
    write(folder/'actual-code-review.json',code);write(folder/'feature-review.json',features)
    write(WORK/'capture-acceptance.json',{'passed':True,'candidate_owner':record(folder/'owner.json'),'candidate_code_review':record(folder/'actual-code-review.json'),
          'candidate_feature_review':record(folder/'feature-review.json'),'stock_code_review':record(stock_folder/'actual-code-review.json'),
          'complete_durable_stock_report_reproduced':True,'fresh_complete55_feature_queue_extension_order_and_chain_recorded':True,
          'fresh_inferior_and_called_JIT_maps_join_all_three_provider_files':True,
          'canonical_stock_raw_feature_arguments_available':False,'canonical_host_Int64_disabled':True,'canonical_unused_default_Float64_enabled':True,
          'original_failed_compiler_wrapper_preserved':True,'successful_restoring_compiler_and_graph_PASS_retained':True,
          'previous_distinct_UDiv_JIT_cost_premise_survived_unproven':True,
          'graph_receipt':graph,'stock_receipt':stock,
          'scope':'One new full-module ordinary loaded/called LLVM/JIT capture, original canary refusal before arithmetic. Complete sections/defined relocations/uncovered bytes checked; undefined callbacks consistency only. No arithmetic/hit-rate/cost/benefit or qualification conclusion; actual cost-premise review remains separate.'})
    print('Candidate ordinary code, actual canonical features and live maps accepted; cost premise unjudged.',flush=True)

if __name__=='__main__':{'pin-runtime':pin_runtime,'review':review}[sys.argv[1]]()
