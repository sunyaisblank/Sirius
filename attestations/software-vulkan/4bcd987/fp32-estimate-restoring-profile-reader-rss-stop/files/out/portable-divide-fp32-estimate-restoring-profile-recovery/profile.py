"""Root-owned stock reuse, full compiler producer, then distinct profile reader.

Successful compiler outputs remain separately sealed if the later reader refuses
an unobserved representation. No provider, science test or benefit is run here.
"""
from pathlib import Path
import hashlib,json,os,re,resource,shutil,stat,subprocess,sys
sys.dont_write_bytecode=True
W=Path(__file__).resolve().parent
R=W.parents[1]
P=json.loads((W/'packet.json').read_text())
T=Path(P['products'])


def validate_action_schema():
    """Required owner contract, checked before creating the operator context."""
    actions=json.loads((W/'actions.json').read_text())
    labels=['reuse-stock','compile-candidate','read-profile']
    assert P['labels']==labels and set(actions)==set(labels)
    required={'kind','revision','argv','frozen_inputs','limits','scope','expected_exit','software_provider'}
    for label in labels:
        action=actions[label]
        assert set(action)==required,(label,set(action))
        assert action['kind']=='build' and action['revision']==P['revision']
        assert type(action['expected_exit']) is int and action['expected_exit']==0
        assert action['software_provider'] is False
        assert isinstance(action['scope'],str) and action['scope']
        assert action['argv']==['/usr/bin/python3','-B',str(W/'profile.py'),'phase',label]
        assert action['limits']==P['limits']=={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
        assert all(type(value) is int and value>0 for value in action['limits'].values())
        files=action['frozen_inputs']
        assert isinstance(files,list) and all(isinstance(name,str) and Path(name).is_absolute() for name in files)
        assert len(files)==len(set(files)) and not set(files)&set(P['outputs'][label])
        assert {str(W/name) for name in ('actions.json','source-freeze.json','operator-context.json')}<=set(files)
    assert (W/'actions.json').read_bytes()==(W/'actions.unpinned.json').read_bytes()
    return actions


def record(path):
    path=Path(path)
    with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
    return {'bytes':path.stat().st_size,'sha256':digest,'mode':stat.S_IMODE(path.stat().st_mode),
            'resolved_path':str(path.resolve(strict=True)),
            'symlink_target':os.readlink(path) if path.is_symlink() else None}


def content(path):return {k:v for k,v in record(path).items() if k in ('bytes','sha256')}


def write(path,value):
    path=Path(path);assert not path.exists(),path
    path.write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')


def static():
    for name,pin in json.loads((W/'source-freeze.json').read_text())['pins'].items():assert record(name)==pin,name


def archived(row):
    archive=Path(row['path']);assert content(archive/'manifest.json')==row['manifest']
    manifest=json.loads((archive/'manifest.json').read_text())['files']
    for member,pin in row['selected_members'].items():
        assert manifest[member]==pin and content(archive/member)==pin,member
        assert not (archive/member).is_symlink()
    return archive


def prior_stop_receipt():
    """Bind the preserved failure without promoting it or replaying staging."""
    spec=P['prior_stop'];archive=archived(spec);old=archive/spec['work_member']
    load=lambda name:json.loads((old/name).read_text())
    owner=load('reuse-stock/owner.json');first=load('root-first-stop.json')
    independent=load('independent-stop-review.json');admission=load('root-owner-stop-admission.json')
    assert owner['status']=='terminal' and owner['passed'] is False and owner['child_exit']==0
    assert owner['source']=={'revision':P['revision'],'status':''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[]
    assert len(owner['observed_births'])==1 and all(row['absent'] for row in owner['observed_births'])
    assert first['passed'] and independent['passed'] and admission['passed']
    assert first['classification']==admission['classification']=='OWNER_METADATA_EXPECTED_EXIT_MISSING'
    assert first['original_owner_passed'] is False and admission['original_owner_passed'] is False
    assert first['original_child_exit']==admission['original_child_exit']==0
    assert first['original_wrapper_exit']==admission['original_wrapper_exit']==1
    assert first['completed_stock_files']==admission['completed_stock_files']==39
    assert not first['candidate_compiler_or_profile_reader_executed'] and not admission['candidate_compiler_or_profile_reader_executed']
    assert len(first['actual_commands'])==2 and [row['exit_code'] for row in first['actual_commands']]==[0,1]
    assert first['actual_commands'][1]['controller_traceback'].endswith("KeyError: 'expected_exit'\n")
    # Archived records retain their historical resolved paths. Only regular
    # members under these two explicit old namespaces are relocated.
    oldwork=Path(spec['original_work']);oldproducts=Path(spec['original_products'])
    def historical(expected,path,original):
        assert not path.is_symlink() and content(path)=={k:expected[k] for k in ('bytes','sha256')}
        assert expected['resolved_path']==str(original)
    for key,name in [('owner','reuse-stock/owner.json'),('output_seal','reuse-stock/output-seal.json'),
                     ('stock_acceptance','stock-reuse-acceptance.json'),('source_freeze','source-freeze.json')]:
        historical(first[key],old/name,oldwork/name)
    assert admission['source_revision']==P['revision']
    historical(admission['root_first_stop'],old/'root-first-stop.json',oldwork/'root-first-stop.json')
    historical(admission['independent_stop_review'],old/'independent-stop-review.json',oldwork/'independent-stop-review.json')
    oldactions=load('actions.json');frozen=load('source-freeze.json')['pins'];before=load('reuse-stock/inputs-before.json')
    assert len(frozen)==316 and len(before)==822
    assert (old/'reuse-stock/inputs-before.json').read_bytes()==(old/'reuse-stock/inputs-after.json').read_bytes()
    assert all('expected_exit' not in action for action in oldactions.values())
    assert set(oldactions)==set(P['labels']) and len(oldactions['reuse-stock']['frozen_inputs'])==319
    for name in oldactions['reuse-stock']['frozen_inputs']:
        path=Path(name);key=str(path.relative_to(R)) if path.is_relative_to(R) else str(path)
        recorded=before[key];assert set(recorded)=={'bytes','sha256','resolved_path'}
        if name in frozen:
            assert recorded=={k:frozen[name][k] for k in recorded}
        else:
            assert path in {oldwork/name for name in ('source-freeze.json','operator-context.json','actions.json')}
            historical(recorded,old/path.relative_to(oldwork),path)
    # Prior root/independent stop audits established all822 current inputs at
    # that stop. Here their exact raw maps/316 metadata are bound, not replayed
    # against retired paths; current recovery source/tool bytes are frozen anew.
    for relative,row in spec['copied_noncontroller_members'].items():
        path=W/relative;member=archive/row['archive_member']
        assert content(path)==content(member)=={k:row[k] for k in ('bytes','sha256')}
    seal=load('reuse-stock/output-seal.json');packet=load('packet.json')
    assert seal['label']=='reuse-stock' and set(seal['outputs'])==set(packet['outputs']['reuse-stock'])
    staged=archive/spec['products_member']/'stock';assert len(list(staged.iterdir()))==39
    copied=load('stock-reuse-acceptance.json');assert copied['passed'] and len(copied['copied39'])==39
    for name,expected in seal['outputs'].items():
        original=Path(name)
        if original.is_relative_to(oldproducts):actual=archive/spec['products_member']/original.relative_to(oldproducts)
        else:
            assert original==oldwork/'stock-reuse-acceptance.json';actual=old/'stock-reuse-acceptance.json'
        historical(expected,actual,original)
    assert copied['copied39']=={str(oldproducts/'stock'/p.name):seal['outputs'][str(oldproducts/'stock'/p.name)] for p in staged.iterdir()}
    assert not (old/'compile-candidate').exists() and not (old/'read-profile').exists()
    return {'manifest':record(archive/'manifest.json'),'original_owner':record(old/'reuse-stock/owner.json'),
            'original_owner_passed':False,'original_child_exit':0,'original_wrapper_exit':1,
            'root_stop':record(old/'root-owner-stop-admission.json'),'independent_stop':record(old/'independent-stop-review.json'),
            'fresh_action_schema_repair_required_and_checked':True,'old_staged_products_adopted_as_new_owner':False,
            'scope':'Immutable causal action-schema stop and completed child0 staging preserved. No success promotion or replay of old822 current paths. Fresh3-owner sequence is separate.'}


def prerequisites():
    prior_stop=prior_stop_receipt()
    arithmetic=archived(P['arithmetic']);aw=arithmetic/P['arithmetic']['work_member']
    result=json.loads((aw/'arithmetic-result.json').read_text())
    terminal=json.loads((aw/'independent-terminal-review.json').read_text())
    root=json.loads((aw/'root-terminal-admission.json').read_text())
    source=json.loads((aw/'root-source-admission.json').read_text())
    assert result['passed'] and terminal['passed'] and root['passed'] and source['passed']
    assert result['original_records']==174268 and result['original_operation_count']==14
    assert result['complete_output_words_per_provider']==697072
    assert result['fresh_candidate_host_and_full_gpu_plus_reused_stock_equal_original_expected']
    assert result['stock_result_reused'] and result['candidate_macro']==P['candidate_macro']
    assert {k:root['arithmetic_result'][k] for k in ('bytes','sha256')}==content(aw/'arithmetic-result.json')
    assert {k:root['independent_terminal_review'][k] for k in ('bytes','sha256')}==content(aw/'independent-terminal-review.json')
    assert {k:result['candidate_header'][k] for k in ('bytes','sha256')}==content(W/'candidate/src/sirius/kernels/portable_binary32.h')
    costarchive=archived(P['cost']);cw=costarchive/P['cost']['work_member']
    cost=json.loads((cw/'cost-review/cost-premise-review.json').read_text())
    admission=json.loads((cw/'root-cost-premise-admission.json').read_text())
    assert cost['passed'] and admission['passed']
    assert cost['decision']==admission['decision']=='COST_PREMISE_SURVIVES_UNPROVEN'
    assert not cost['retained_arithmetic_executed_by_canary'] and cost['no_arithmetic_or_benefit_admission']
    assert {k:admission['cost_review'][k] for k in ('bytes','sha256')}==content(cw/'cost-review/cost-premise-review.json')
    from source_profile import kernel_source_profile
    source_profile=kernel_source_profile(P,R,W,record,content)
    return {'prior_stop':prior_stop,'arithmetic_result':record(aw/'arithmetic-result.json'),'arithmetic_terminal':record(aw/'independent-terminal-review.json'),
            'arithmetic_root_admission':record(aw/'root-terminal-admission.json'),'arithmetic_manifest':record(arithmetic/'manifest.json'),
            'cost':record(cw/'cost-review/cost-premise-review.json'),'cost_root_admission':record(cw/'root-cost-premise-admission.json'),
            'source_profile':source_profile,
            'scope':'Exact pushed accepted finite prerequisite receipts and relevant source/profile reuse. No original numeric/cost workload repeated and no production adoption.'}


def terminal(folder):
    folder=Path(folder);owner=json.loads((folder/'owner.json').read_text())
    assert owner['status']=='terminal' and owner['passed'] and owner['child_exit']==0
    assert owner['stop_reason'] is None and owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[]
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert (folder/'inputs-before.json').read_bytes()==(folder/'inputs-after.json').read_bytes()
    return owner


def context():
    validate_action_schema()
    saved=json.loads((W/'operator-context.json').read_text())
    for key,name in [('source_freeze','source-freeze.json'),('actions','actions.json'),('preparation','preparation.json')]:assert saved[key]==record(W/name)
    assert saved['revision']==P['revision']
    assert (W/'actions.json').read_bytes()==(W/'actions.unpinned.json').read_bytes()
    return saved


def seal(label):
    write(W/label/'output-seal.json',{'label':label,'recipe':P['recipes'][label],
          'source_freeze':record(W/'source-freeze.json'),'operator_context':record(W/'operator-context.json'),
          'inputs_before':record(W/label/'inputs-before.json'),
          'outputs':{name:record(name) for name in P['outputs'][label]}})


def accepted(label):
    owner=terminal(W/label);action=json.loads((W/'actions.json').read_text())[label]
    assert owner['source']=={'revision':P['revision'],'status':''} and owner['argv']==action['argv'] and owner['limits']==P['limits']
    before=json.loads((W/label/'inputs-before.json').read_text())
    for name in action['frozen_inputs']:
        path=Path(name);key=str(path.relative_to(R)) if path.is_relative_to(R) else str(path)
        # Mature owner identity has exactly the original three fields. Rich
        # mode/symlink metadata remains separately bound by static().
        assert set(before[key])=={'bytes','sha256','resolved_path'}
        assert before[key]=={k:record(path)[k] for k in before[key]},name
    saved=json.loads((W/label/'output-seal.json').read_text())
    assert saved['label']==label and saved['recipe']==P['recipes'][label]
    for key,path in [('source_freeze',W/'source-freeze.json'),('operator_context',W/'operator-context.json'),('inputs_before',W/label/'inputs-before.json')]:assert saved[key]==record(path)
    assert saved['outputs']=={name:record(name) for name in P['outputs'][label]}
    return {'owner':record(W/label/'owner.json'),'inputs_before':record(W/label/'inputs-before.json'),
            'inputs_after':record(W/label/'inputs-after.json'),'output_seal':record(W/label/'output-seal.json')}


def stock_receipt():
    archive=archived(P['stock_archive']);work=Path(P['historical_stock_producer']['work'])
    sourcework=Path(P['original_stock_work'])
    assert content(work/'historical-raw-manifest.json')==P['historical_stock_producer']['raw_archive_manifest']
    historical_manifest=json.loads((work/'historical-raw-manifest.json').read_text())['files']
    oldroot=Path(P['historical_stock_producer']['original_work_path'])
    for relative in ('compile-products-stock/owner.json','compile-products-stock/inputs-before.json',
                     'compile-products-stock/inputs-after.json','products-acceptance.json'):
        member='files/'+str((oldroot/relative).relative_to(R))
        assert content(work/relative)==historical_manifest[member]
    owner=terminal(work/'compile-products-stock')
    before=json.loads((work/'compile-products-stock/inputs-before.json').read_text())
    assert owner['argv']==P['original_stock_recipe']
    accepted_products=json.loads((work/'products-acceptance.json').read_text());assert accepted_products['passed']
    assert accepted_products['original_builder_owners']['stock']==content(work/'compile-products-stock/owner.json')
    for key in P['stock_reuse_tools']:
        tool=Path(key);old=[row for row in before.values() if row['resolved_path']==str(tool.resolve(strict=True))]
        assert old and all(content(tool)=={k:v[k] for k in ('bytes','sha256')} for v in old),key
    original=sourcework/'stock/src/sirius/kernels'
    for file in original.iterdir():
        key=str((R/'out/computed-carry-exit-original-prerequisites/stock/src/sirius/kernels'/file.name).relative_to(R))
        assert content(file)=={k:before[key][k] for k in ('bytes','sha256')}
    builderkey='out/computed-carry-exit-original-prerequisites/stock/scripts/build-retained-kernels.py'
    assert content(sourcework/'stock/scripts/build-retained-kernels.py')=={k:before[builderkey][k] for k in ('bytes','sha256')}
    products=Path(P['stock_reuse_directory']);canonical=R/'bin/linux-gcc/src/sirius/backend/retained'
    expected={row['module']:row['stock'] for row in accepted_products['modules']}
    assert len(expected)==35
    for path in products.iterdir():
        assert path.read_bytes()==(canonical/path.name).read_bytes(),path
        if path.suffix=='.spv':assert content(path)==expected[path.name]
    assert len(list(products.iterdir()))==39
    return {'manifest':record(archive/'manifest.json'),'original_owner':record(work/'compile-products-stock/owner.json'),
            'original_before':record(work/'compile-products-stock/inputs-before.json'),
            'original_after':record(work/'compile-products-stock/inputs-after.json'),
            'original_products_acceptance':record(work/'products-acceptance.json'),
            'raw_historical_receipt_scope':P['historical_stock_producer'],
            'recipe_source33_tools5_and_canonical39_equal':True,
            'scope':'Relevant byte-producing original33 source/builder/five tool/recipe identities plus39 products and terminal input-equality. No claim every unrelated historical owner input matches current; exact14 consumed-source closure is separately checked.'}


def saved_assembly_receipt():
    archive=archived(P['assembly_archive']);work=archive/P['assembly_archive']['work_member']
    result=json.loads((work/'profile-result.json').read_text())
    acceptance=json.loads((work/'profile-acceptance.json').read_text())
    independent=json.loads((work/'independent-terminal-review.json').read_text())
    assert result['passed'] and acceptance['passed'] and independent['passed']
    assert result['original_failed_wrapper_passed'] is False and result['original_failed_wrapper_child_exit']==1
    assert {k:acceptance['profile_result'][k] for k in ('bytes','sha256')}==content(work/'profile-result.json')
    rows={row['module']:row for row in result['profile']['modules']};assert len(rows)==35
    for name in P['module_arrays']:
        row=rows[name]['assembly_binary_joins']['stock'];saved=P['saved_assembly']['stock/'+name]
        assert content(saved['path'])==saved['pin']=={k:row['assembly'][k] for k in ('bytes','sha256')}
        assert content(Path(P['stock_reuse_directory'])/name)=={k:row['source_module'][k] for k in ('bytes','sha256')}
        # Historical resolved identities are preserved, with only this explicit
        # archive namespace relocation for the32 fresh stock disassemblies.
        old=Path(row['assembly']['resolved_path']);actual=Path(saved['path'])
        assert actual==old or (old.is_relative_to(R/'out/portable-divide-fp64-estimate-profile-recovery') and
            actual==work/old.relative_to(R/'out/portable-divide-fp64-estimate-profile-recovery'))
    return {'accepted_original_profile':record(work/'profile-acceptance.json'),
            'complete35_saved_stock_assembly_binary_pairs_equal':True,
            'reuse_count_with18_unchanged_candidate_aliases':53,'new_candidate_disassemblies':17,
            'old_failed_profile_owner_remains_false_child1':True}


def pin_context():
    validate_action_schema();static();prerequisites();stock_receipt();saved_assembly_receipt()
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()==P['revision']
    assert not subprocess.check_output(['git','status','--porcelain=v1'],cwd=R)
    assert not T.exists() and not (W/'analysis').exists() and not (W/'operator-context.json').exists()
    assert (W/'actions.json').read_bytes()==(W/'actions.unpinned.json').read_bytes()
    for label in P['labels']:assert not (W/label).exists()
    for files in P['outputs'].values():assert all(not Path(name).exists() for name in files)
    write(W/'operator-context.json',{'revision':P['revision'],'source_freeze':record(W/'source-freeze.json'),
          'actions':record(W/'actions.json'),'preparation':record(W/'preparation.json')})


def reuse_stock():
    proof=stock_receipt();target=T/'stock';assert not target.exists();target.mkdir(parents=True)
    for path in Path(P['stock_reuse_directory']).iterdir():shutil.copyfile(path,target/path.name)
    write(W/'stock-reuse-acceptance.json',{'passed':True,'historical':proof,'copied39':{str(p):record(p) for p in target.iterdir()}})


def compile_candidate():
    accepted('reuse-stock');assert not (T/'candidate').exists()
    for argv in P['recipes']['compile-candidate']:subprocess.run(argv,cwd=R,check=True)
    candidate=T/'candidate'
    assert {str(p) for p in candidate.iterdir()}==set(P['candidate_products'])
    write(W/'compiler-products.json',{'schema':'sirius-fp32-restoring-full-profile-compiler-products-v1','passed':True,
          'products':{name:record(name) for name in P['candidate_products']},'selected17':P['changed_modules'],
          'compiler_policy_and_validation_completed':True,'structural_graph_accepted':False,
          'scope':'Complete original35 candidate recipes; selected17 raw/assembler products retained. Separate profile reader may refuse without changing this compiler result.'})


def read_profile():
    accepted('reuse-stock');accepted('compile-candidate')
    assert json.loads((W/'compiler-producer-acceptance.json').read_text())['passed']
    saved=saved_assembly_receipt()
    from products_profile import products_review
    result=products_review(P,R,W,T,W)
    assert result['changed_modules']==17 and result['unchanged_modules']==18
    assert len(list((W/'analysis/candidate').glob('*.spvasm')))==17
    write(W/'profile-result.json',{'schema':'sirius-fp32-restoring-complete-retained-profile-v1','passed':True,
          'profile':result,'saved_assembly_reuse':saved,'prerequisites':prerequisites(),'arithmetic_or_benefit_or_adoption_admitted':False})


def phase(label):
    validate_action_schema()
    resource.setrlimit(resource.RLIMIT_CORE,(0,0));_,prior=resource.getrlimit(resource.RLIMIT_FSIZE)
    cap=64*1024**2;limit=cap if prior==resource.RLIM_INFINITY else min(prior,cap);resource.setrlimit(resource.RLIMIT_FSIZE,(limit,limit))
    static();context();prerequisites()
    for name in P['outputs'][label]:assert not Path(name).exists(),name
    {'reuse-stock':reuse_stock,'compile-candidate':compile_candidate,'read-profile':read_profile}[label]()
    static();seal(label)


def review(label):
    validate_action_schema()
    static();context();prerequisites();receipt=accepted(label)
    if label=='compile-candidate':
        write(W/'compiler-producer-acceptance.json',{'schema':'sirius-fp32-restoring-full-profile-compiler-acceptance-v1',
              'passed':True,'compiler_result':record(W/'compiler-products.json'),'compiler_receipt':receipt,
              'structural_graph_accepted':False,'provider_or_numerical_or_benefit_executed':False})
    elif label=='read-profile':
        result=json.loads((W/'profile-result.json').read_text());assert result['passed']
        assert result['prerequisites']==prerequisites()
        write(W/'profile-acceptance.json',{'schema':'sirius-fp32-restoring-complete-retained-profile-acceptance-v1',
              'passed':True,'profile_result':record(W/'profile-result.json'),'reader_receipt':receipt,
              'compiler_acceptance':record(W/'compiler-producer-acceptance.json'),'compiler_receipt':accepted('compile-candidate'),
              'stock_receipt':accepted('reuse-stock'),'all35_modules_40main_2owner_arrays_7tables_accepted':True,
              'named_graphs':13,'inlined_camera_footprints':4,'arithmetic_or_integrated_or_benefit_or_adoption_admitted':False,
              'camera_scope':'Guarded restoring dataflow and exact integer/loop/interface inventories; full normalization/Round/fallback equivalence additionally relies on frozen source and original O1 inlining.'})


if __name__=='__main__':
    if sys.argv[1]=='pin-context':pin_context()
    elif sys.argv[1]=='phase':phase(sys.argv[2])
    elif sys.argv[1]=='review':review(sys.argv[2])
    else:raise AssertionError('unknown phase')
