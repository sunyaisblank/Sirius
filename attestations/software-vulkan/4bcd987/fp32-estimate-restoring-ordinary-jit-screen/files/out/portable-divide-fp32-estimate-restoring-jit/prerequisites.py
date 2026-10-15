"""Archived restoring graph/compiler truth and original canonical stock producer receipts."""
from pathlib import Path
import hashlib,importlib.util,json,os,re,shlex,stat

WORK=Path(__file__).resolve().parent
ROOT=WORK.parents[1]
P=json.loads((WORK/'inputs.json').read_text())

def load(path):return json.loads(Path(path).read_text())
def record(path):
    path=Path(path)
    with path.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
    return {'bytes':path.stat().st_size,'sha256':digest,'resolved_path':str(path.resolve(strict=True)),
            'mode':stat.S_IMODE(path.stat().st_mode),'symlink_target':os.readlink(path) if path.is_symlink() else None}
def content(path):
    value=record(path);return {k:value[k] for k in ('bytes','sha256')}
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result

def accepted(folder):
    owner=load(folder/'owner.json')
    assert owner['status']=='terminal' and owner['passed'] is True and owner['child_exit']==0
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and owner['stop_reason'] is None
    assert owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[] and owner['observed_births'] and all(b['absent'] for b in owner['observed_births'])
    assert (folder/'inputs-before.json').read_bytes()==(folder/'inputs-after.json').read_bytes()
    return owner

def original_path(name):
    path=Path(name);return path if path.is_absolute() else ROOT/path

def graph_joined(name,pin):
    original=original_path(name);retired=ROOT/P['graph_original_out']
    actual=Path(P['graph_archive'])/'files'/original.relative_to(ROOT) if original.is_relative_to(retired) else original
    seen=record(actual)
    if actual!=original:
        assert not actual.is_symlink() and pin['resolved_path']==str(original),name
        seen['resolved_path']=str(original)
    assert {k:seen[k] for k in pin}==pin,name
    return actual

def verify_graph():
    assert re.fullmatch(r'[0-9a-f]{40}',P['graph_commit']), 'actual immutable restoring graph required'
    archive=Path(P['graph_archive']);manifest=archive/'manifest.json'
    assert content(manifest)==P['graph_manifest']
    rows=load(manifest)['files']
    for name,pin in rows.items():assert not (archive/name).is_symlink() and content(archive/name)==pin,name
    old=archive/'files'/P['graph_original_out'];packet=load(old/'packet.json')
    assert packet['compiler_commit']==P['compiler_commit'] and packet['compiler_manifest']==P['compiler_manifest']
    assert packet['compiler_archive']==P['compiler_archive'] and packet['candidate_header_sha256']==P['candidate_header_sha256']
    # This pure completed-compiler verifier joins the current102/636 producer,
    # not prior JIT933, graph/control algorithms, tools or provider execution.
    prior=module('immutable_restoring_compiler_prerequisite_reader',old/'proof.py');prior.R=ROOT
    chain=prior.completed()
    for name,pin in load(old/'source-freeze.json')['pins'].items():graph_joined(name,pin)
    context=load(old/'operator-context.json');assert context['revision']==P['operator_revision']
    for key,name in [('actions','actions.json'),('source_freeze','source-freeze.json'),('preparation','preparation.json')]:graph_joined(ROOT/P['graph_original_out']/name,context[key])
    owner=accepted(old/'read-graph');assert owner['source']=={'revision':P['operator_revision'],'status':''}
    action=load(old/'actions.json')['read-graph'];assert owner['argv']==action['argv'] and owner['limits']==P['limits']==action['limits']
    before=load(old/'read-graph/inputs-before.json')
    for name,pin in before.items():graph_joined(name,pin)
    seal=load(old/'read-graph/output-seal.json')
    for key,name in [('source_freeze','source-freeze.json'),('operator_context','operator-context.json'),('inputs_before','read-graph/inputs-before.json')]:graph_joined(ROOT/P['graph_original_out']/name,seal[key])
    expected={str(ROOT/P['graph_original_out']/name):pin for name,pin in ((Path(key).name,value) for key,value in seal['outputs'].items())}
    assert set(expected)=={str(ROOT/P['graph_original_out']/name) for name in ('graph-result.json','parser-control-result.json')}
    for name,pin in expected.items():graph_joined(name,pin)
    result=load(old/'graph-result.json');tests=load(old/'parser-control-result.json');approval=load(old/'graph-acceptance.json')
    assert result['schema']=='sirius-fp32-estimate-restoring-full-endpoint-graph-result-v1'
    assert approval['schema']=='sirius-fp32-estimate-restoring-full-endpoint-graph-acceptance-v1'
    assert result['passed'] and tests['passed'] and approval['passed'] and approval['structural_graph_accepted']
    assert tests['positive_count']==3 and tests['negative_count']==packet['required_graph_negative_controls']
    assert tests['helper_normalization_controls']['passed'] and tests['helper_normalization_controls']['negative_count']==2
    assert not result['arithmetic_execution'] and not result['cost_or_hit_rate_or_adoption_established'] and not approval['arithmetic_or_cost_or_adoption_admitted']
    for key,name in [('graph_owner','read-graph/owner.json'),('output_seal','read-graph/output-seal.json'),('graph_result','graph-result.json'),('parser_controls','parser-control-result.json')]:graph_joined(ROOT/P['graph_original_out']/name,approval[key])
    for saved in (result['completed_compiler_chain'],approval['immutable_compiler_chain']):
        assert {k:v for k,v in saved.items() if k!='fresh_saved_birth_absence'}=={k:v for k,v in chain.items() if k!='fresh_saved_birth_absence'}
        assert all(b['saved_birth_absent'] for values in saved['fresh_saved_birth_absence'].values() for b in values)
    assert approval['original_failed_passed'] is False and approval['original_failed_child_exit']==1
    assert result['candidate_macro']=='SIRIUS_PB32_DIVIDE_FP32_ESTIMATE=1' and result['candidate_header_sha256']==P['candidate_header_sha256']
    assert result['unchanged_normalized_PB32_helper_instances']['unchanged_non_Divide_helper_bodies']==23
    assert result['assembler_and_readback_complete_instructions_equal']==packet['complete_module_instruction_count']
    actual=result['actual_Divide']
    for flag in ('original_complete_prefix_normalized_equal','original_complete_fallback_suffix_normalized_equal','exact_new_region_instruction_and_block_coverage','NoUnsignedWrap_and_NoSignedWrap_absent','corrective_UDiv_absent','only_original_fallback_UDiv_remains'):assert actual[flag]
    assert [step['shift'] for step in actual['restoring_steps']]==[4,3,2,1,0]
    for name in ('independent-terminal-review.json','root-terminal-admission.json'):
        assert load(old/name)['passed'],name
    candidate=Path(chain['candidate_products']['modules/candidate-endpoint-normal.spv']['resolved_path'])
    assert content(candidate)==P['candidate_module']
    return {'graph_commit':P['graph_commit'],'graph_manifest':record(manifest),'graph_acceptance':record(old/'graph-acceptance.json'),
            'graph_result':record(old/'graph-result.json'),'graph_owner':record(old/'read-graph/owner.json'),'graph_output_seal':record(old/'read-graph/output-seal.json'),
            'graph_terminal_review':record(old/'independent-terminal-review.json'),'candidate_module':record(candidate),
            'original_failed_owner':chain['original_failed_owner'],'original_failed_passed':False,'original_failed_child_exit':1,
            'compiler_owner':chain['compiler_owner'],'compiler_products':chain['compiler_products'],'compiler_acceptance':chain['compiler_acceptance'],
            'graph_passed':True,'restoring_compiler_passed':True,'graph_or_policy_or_compiler_repeated':False,
            'previous_distinct_UDiv_JIT_cost_premise_survived_unproven':True,
            'fresh_graph_birth_absence':prior.fresh_births(owner),'arithmetic_execution':False}

def verify_stock():
    data=load(WORK/'stock-reuse-inputs.json');archive=Path(data['archive']);history=Path(data['historical_work'])
    manifest_path=archive/data['manifest_name'];assert content(manifest_path)['sha256']==data['manifest_sha256']
    current=load(manifest_path)['files'];original_manifest=Path(data['original_manifest'])
    assert content(original_manifest)['sha256']==data['original_manifest_sha256']
    original=load(original_manifest)
    assert original['source_revision']==P['product_revision'] and original['evidence_commit']==data['original_evidence_commit'] and original['evidence_tag']==data['original_evidence_tag']
    used={}
    def manifested(path):
        path=Path(path);member=str(path.relative_to(archive));assert not path.is_symlink() and content(path)==current[member],member
        prefix='stock-comparison/files/'
        if member.startswith(prefix):
            old=member[len(prefix):];assert content(path)==original['files'][old],old
            used[old]=content(path)
    manifested(original_manifest)
    recovery=load(history/'input-recovery.json');manifested(history/'input-recovery.json')
    saved={(row['input'],row['sha256']):row for row in recovery['records']}
    lookup={(row['original_key'],row['pin']['sha256']):row for row in data['joins']};assert len(lookup)==len(data['joins'])==638
    consumed=set()
    def joined(key,pin,phase):
        row=lookup[(key,pin['sha256'])];assert row['pin']==pin and phase in row['phases']
        original=original_path(key);actual=Path(row['actual']);assert row['original_path']==str(original)
        assert content(actual)=={k:pin[k] for k in ('bytes','sha256')}
        method=row['method']
        if method=='exact_current_identity':assert actual==original and str(actual.resolve(strict=True))==pin['resolved_path']
        elif method=='whole_byte_manifest_archived_owned_file':
            assert actual==archive/'stock-comparison/files'/row['manifest_member'] and original==ROOT/row['manifest_member']
            assert not actual.is_symlink() and pin['resolved_path']==str(original);manifested(actual)
        elif method=='exact_original_input_recovery':
            old=saved[(key,pin['sha256'])];assert {k:old[k] for k in pin}==pin
            assert actual==archive/'stock-comparison/files'/old['saved_path'] and row['manifest_member']==old['saved_path'];manifested(actual)
        elif method=='selected_canonical_producer':
            assert row['original_method']=='exact_current_identity' and pin['resolved_path']==str(original)
            assert key in ('bin/linux-gcc/src/sirius/backend/CMakeFiles/sirius_backend.dir/vulkan/vulkan_device.cpp.o','bin/linux-gcc/src/sirius/libsirius_base.a')
            assert actual==archive/row['selected_manifest_member'];manifested(actual)
        elif method=='exact_historical_Git_blob':
            assert row['git_commit']==P['product_revision'] and row['git_path']==str(original.relative_to(ROOT))
            raw=actual.read_bytes();assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==row['git_blob']
            assert actual==WORK/'recovered-source'/row['git_path'] and pin['resolved_path']==str(original)
        else:raise AssertionError(method)
        consumed.add((key,pin['sha256']));return actual
    owners={}
    for phase,count in (('compile-probe',545),('link-probe',547),('capture-stock-endpoint-normal',633)):
        folder=history/phase;owners[phase]=accepted(folder)
        assert owners[phase]['source']=={'revision':P['product_revision'],'status':''} and owners[phase]['limits']==P['limits']
        for name in ('owner.json','inputs-before.json','inputs-after.json'):manifested(folder/name)
        before=load(folder/'inputs-before.json');assert len(before)==count
        for key,pin in before.items():joined(key,pin,phase)
    assert consumed==set(lookup)
    compile_owner=owners['compile-probe'];link=owners['link-probe']
    source=str(ROOT/'out/guarded-product-high-runtime/short_input_jit_probe.cpp')
    assert source in compile_owner['argv'] and '-O3' in compile_owner['argv'] and '-MMD' in compile_owner['argv']
    device='bin/linux-gcc/src/sirius/backend/CMakeFiles/sirius_backend.dir/vulkan/vulkan_device.cpp.o';base='bin/linux-gcc/src/sirius/libsirius_base.a'
    assert link['argv'].count(str(ROOT/device))==link['argv'].count(str(ROOT/base))==1 and not any('libsirius_backend.a' in arg for arg in link['argv'])
    stock_host=archive/'stock-comparison/files/bin/linux-gcc/guarded-product-high-runtime/short_input_jit_probe'
    stock_module=archive/'stock-comparison/files/bin/linux-gcc/guarded-product-high-runtime/modules/stock-endpoint-normal.spv'
    canary=Path(P['canary']);probe=history/'short_input_jit_probe.cpp'
    for path in (stock_host,stock_module,canary,probe,history/'preparation.json'):manifested(path)
    assert content(canary)==content(stock_host)==P['canary_identity'] and content(stock_module)==P['stock_module']
    assert content(probe)=={'bytes':4212,'sha256':'853c361fa68e8afd3c532cf7cf5847837f1620cc047bd773a39d0955695241d7'}
    link_before=load(history/'link-probe/inputs-before.json')
    dev=joined(device,link_before[device],'link-probe');lib=joined(base,link_before[base],'link-probe')
    assert content(dev)=={'bytes':291264,'sha256':'181071b73bc05e54ea60116f654aa6f45c61451a0e422357c638ed1b105131be'}
    assert content(lib)=={'bytes':35862,'sha256':'1fb43165640e4059b46e8cd2802a954e529788cdc17e832fa1c583f920014d76'}
    assert content(ROOT/'src/sirius/backend/vulkan/vulkan_device.cpp')=={'bytes':56500,'sha256':'c3765cbb0e34bcd509ea63da6f264fbe4e34d2a6b4584bc0bdeedf1a75cf77ae'}
    assert content(ROOT/'src/sirius/backend/device.h')=={'bytes':9669,'sha256':'be88c5cbec6c8dadabf10a9feb52e4df562ad811fcdd41aa90b4645852616983'}
    stock=history/'capture-stock-endpoint-normal'
    for path in sorted(stock.iterdir()):
        if path.is_file():manifested(path)
    assert len(used)==56 and len([name for name in used if '/capture-stock-endpoint-normal/' in name])==33
    report=load(stock/'actual-code-review.json')
    assert report['actual_loaded_module']['sha256']==P['stock_module']['sha256'] and report['single_completed_dispatch']['completed']
    assert report['ordinary_pipelines']==P['ordinary_pipelines'] and report['ordinary_ir_export']['sha256']=='3271b6d61996590a4c74f822659bbe8c33cd3cb0b91f4bdf836ae71a95a290b7'
    return {'canonical_stock_host':record(canary),'stock_module':record(stock_module),'stock_source':record(probe),'canonical_device_object':record(dev),'canonical_base_archive':record(lib),
            'original_compile_owner':record(history/'compile-probe/owner.json'),'original_link_owner':record(history/'link-probe/owner.json'),
            'original_stock_capture_owner':record(stock/'owner.json'),'stock_code_review':record(stock/'actual-code-review.json'),
            'selected_archive_commit':data['evidence_commit'],'selected_archive_manifest':record(manifest_path),
            'whole_byte_original_manifest':record(original_manifest),'input_recovery':record(history/'input-recovery.json'),'matched_historical_input_versions':638,
            'selected_manifest_members':len(used),'selected_stock_raw_members':33,'stock_capture_folder':str(stock),
            'host_is_unmodified_canonical_Create':True,'host_rebuilt':False,'stock_recaptured':False,
            'historical_stock_raw_feature_argument_capture_available':False,
            'scope':'Direct638 original input-version/13Git-source/56selectedmember/33raw stock and canonical .o/base byte joins to immutable0ed; not a replay of prior933 and not a whole older multicapture backup claim. Original canonical unusedFloat64/Int64zero host retained.'}
