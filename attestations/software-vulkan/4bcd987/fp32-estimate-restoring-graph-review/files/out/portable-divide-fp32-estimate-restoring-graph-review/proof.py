"""Immutable successful compiler receipts; no compiler/policy/provider operation."""
from pathlib import Path
import hashlib,json,os,re,stat

W=Path(__file__).resolve().parent
R=W.parents[1]
P=json.loads((W/'packet.json').read_text())

def record(path):
    path=Path(path)
    with path.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
    return {'bytes':path.stat().st_size,'sha256':h,'mode':stat.S_IMODE(path.stat().st_mode),
            'resolved_path':str(path.resolve(strict=True)),'symlink_target':os.readlink(path) if path.is_symlink() else None}

def content(path):
    value=record(path);return {k:value[k] for k in ('bytes','sha256')}

def archive():return Path(P['compiler_archive'])
def oldwork():return archive()/'files'/P['compiler_out']
def products():return archive()/'files'/P['compiler_bin']
def original_path(name):
    path=Path(name);return R/path if not path.is_absolute() else path
def relocated(name):
    path=original_path(name)
    return archive()/'files'/path.relative_to(R) if any(path.is_relative_to(R/P[key]) for key in ('compiler_out','compiler_bin')) else path
def joined(name,pin):
    original=original_path(name);actual=relocated(name);seen=record(actual)
    if actual!=original:
        assert not actual.is_symlink() and pin['resolved_path']==str(original),name
        seen['resolved_path']=str(original)
    assert {k:seen[k] for k in pin}==pin,name
    return actual
def historical_record(name):
    original=original_path(name);seen=record(relocated(name))
    if relocated(name)!=original:seen['resolved_path']=str(original)
    return seen
def fresh_births(owner):
    boot=Path('/proc/sys/kernel/random/boot_id').read_text().strip();rows=[]
    for birth in owner['observed_births']:
        assert birth['absent'];path=Path('/proc')/str(birth['pid'])/'stat'
        try:
            text=path.read_text();fields=text[text.rfind(')')+2:].split();ticks=fields[19]
        except (FileNotFoundError,ProcessLookupError):ticks=None
        assert not (boot==owner['boot_id'] and ticks==birth['saved_start_ticks']),birth
        rows.append({'pid':birth['pid'],'saved_start_ticks':birth['saved_start_ticks'],'current_start_ticks':ticks,
                     'same_boot':boot==owner['boot_id'],'saved_birth_absent':True})
    return rows

def completed():
    assert re.fullmatch(r'[0-9a-f]{40}',P['compiler_commit']), 'immutable compiler archive required'
    a=archive();manifest=a/'manifest.json';assert content(manifest)==P['compiler_manifest']
    members=json.loads(manifest.read_text())['files'];assert len(members)==P['compiler_manifest_members']
    for name,pin in members.items():
        path=a/name;assert not path.is_symlink() and content(path)==pin,name
    old=oldwork();frozen=json.loads((old/'source-freeze.json').read_text())
    assert sum(len(rows) for rows in frozen.values())==P['compiler_static_count']==102
    for name,pin in frozen['private_inputs'].items():joined(R/P['compiler_out']/name,pin)
    for name,pin in frozen['external_inputs'].items():joined(name,pin)
    context=json.loads((old/'operator-context.json').read_text());assert context['revision']==P['operator_revision'] and context['status']==''
    for key,name in [('actions','actions.json'),('source_freeze','source-freeze.json'),('preparation','preparation.json')]:joined(R/P['compiler_out']/name,context[key])
    actions=json.loads((old/'actions.json').read_text());unpinned=json.loads((old/'actions.unpinned.json').read_text())
    assert list(actions)==list(unpinned)==['compile-endpoint-normal']
    action=actions['compile-endpoint-normal'];assert action['revision']==P['operator_revision']
    assert unpinned['compile-endpoint-normal']['revision']=='PENDING_OPERATOR_CLEAN_CONTEXT'
    unpinned['compile-endpoint-normal']['revision']=P['operator_revision'];assert unpinned==actions
    folder=old/'compile-endpoint-normal';owner=json.loads((folder/'owner.json').read_text())
    assert owner['status']=='terminal' and owner['passed'] is True and owner['child_exit']==0
    assert owner['stop_reason'] is None and owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[]
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['source']=={'revision':P['operator_revision'],'status':''} and owner['limits']==P['limits'] and owner['argv']==action['argv']
    assert owner['observed_births'] and all(b['absent'] for b in owner['observed_births'])
    assert (folder/'inputs-before.json').read_bytes()==(folder/'inputs-after.json').read_bytes()
    before=json.loads((folder/'inputs-before.json').read_text());assert len(before)==P['compiler_consumed_count']==636
    for name,pin in before.items():joined(name,pin)
    seal=json.loads((folder/'output-seal.json').read_text());result=json.loads((old/'compiler-products.json').read_text())
    for key,name in [('source_freeze','source-freeze.json'),('operator_context','operator-context.json'),('inputs_before','compile-endpoint-normal/inputs-before.json'),('recipe','compiler-recipes.json')]:joined(R/P['compiler_out']/name,seal[key])
    expected={str(R/P['compiler_bin']/n):historical_record(R/P['compiler_bin']/n) for n in P['product_names']}
    expected[str(R/P['compiler_out']/'compiler-products.json')]=historical_record(R/P['compiler_out']/'compiler-products.json')
    assert seal['label']=='compile-endpoint-normal' and seal['outputs']==expected and seal['structural_graph_accepted'] is False
    assert result['passed'] and result['original_compile_ssa_assemble_validate_complete']
    assert result['stock_equals_current_canonical'] and result['private_single_FDiv32_policy_and_original_normal_policy_passed']
    assert result['exact_one_NoContraction_and_original_two_RTE_lines_inverse_equal_raw']
    assert result['complete_descriptor_and_workgroup_type_decorations_equal'] and result['execution_modes_equal'] and result['ordered_barrier_operands_equal']
    assert not result['structural_graph_accepted'] and not result['runtime_or_numerical_admission'] and not result['benefit_established']
    for key,name in [('raw_post_SSA_module',P['product_names'][0]),('raw_preannotation_assembly',P['product_names'][1]),('actual_assembler_input',P['product_names'][2])]:joined(R/P['compiler_bin']/name,result[key])
    for key,name in [('module',P['product_names'][3]),('assembly',P['product_names'][4])]:joined(R/P['compiler_bin']/name,result['modules']['candidate'][key])
    assert result['modules']['candidate']['module']['sha256']==P['candidate_module_sha256']
    header=old/'candidate/src/sirius/kernels/portable_binary32.h';assert content(header)['sha256']==P['candidate_header_sha256']
    # Reuse actual stock/source/tool/guard receipts without re-executing old policy.
    reuse=json.loads((old/'stock-reuse.json').read_text());stock_archive=Path(reuse['archive'])
    assert content(stock_archive/'manifest.json')['sha256']==reuse['manifest_sha256']
    sm=json.loads((stock_archive/'manifest.json').read_text())['files']
    for name,pin in reuse['archive_pins'].items():
        assert content(name)==pin,name
        if Path(name)!=stock_archive/'manifest.json':assert sm[str(Path(name).relative_to(stock_archive))]==pin
    history=Path(reuse['historical_work']);prior=json.loads((history/'compile-endpoint-normal/owner.json').read_text())
    assert prior['passed'] and prior['status']=='terminal' and prior['child_exit']==0 and prior['stop_reason'] is None
    assert prior['source_unchanged'] and prior['whole_inputs_unchanged'] and prior['cleanup_errors']==[] and prior['remaining_owned_processes']==[]
    assert (history/'compile-endpoint-normal/inputs-before.json').read_bytes()==(history/'compile-endpoint-normal/inputs-after.json').read_bytes()
    sb=json.loads((history/'compile-endpoint-normal/inputs-before.json').read_text())
    assert len(reuse['source_joins'])==8 and len(reuse['tool_joins'])==16
    for row in reuse['source_joins']:
        pin={k:sb[row['original_input_key']][k] for k in ('bytes','sha256')}
        assert all(content(relocated(row[key]))==pin for key in ('historical','private','canonical')),row
    for row in reuse['tool_joins']:
        pin=sb[row['original_input_key']];assert content(row['path'])=={k:pin[k] for k in ('bytes','sha256')}
        assert str(Path(row['path']).resolve(strict=True))==pin['resolved_path']
    assert Path(reuse['stock_module']).read_bytes()==(R/'bin/linux-gcc/src/sirius/backend/retained/retained_endpoint_portable_normal_sum.spv').read_bytes()
    assert result['stock_reuse']['source_joins']==8 and result['stock_reuse']['tool_joins']==16
    for key,name in [('module','stock_module'),('assembly','stock_assembly')]:
        assert result['modules']['stock'][key]==record(reuse[name])
        joined(result['modules']['stock'][key]['resolved_path'],result['modules']['stock'][key])
    for key,name in [('original_owner','compile-endpoint-normal/owner.json'),('original_acceptance','compiler-acceptance.json')]:assert result['stock_reuse'][key]==record(history/name)
    sa=json.loads((history/'compiler-acceptance.json').read_text());sr=json.loads((history/'compiler-result.json').read_text())
    assert sa['passed'] and sa['compiler_owner']['sha256']==content(history/'compile-endpoint-normal/owner.json')['sha256']
    assert sa['compiler_result']['sha256']==content(history/'compiler-result.json')['sha256']
    assert sr['stock_equals_current_canonical'] and sr['stock_builder_policy_unchanged']
    assert result['stock_reuse']['archive_manifest']==record(stock_archive/'manifest.json') and result['stock_reuse']['guard_controls']==record(reuse['guard_controls'])
    acceptance=json.loads((old/'compiler-producer-acceptance.json').read_text());assert acceptance['passed'] and acceptance['schema']=='sirius-fp32-estimate-compiler-producer-acceptance-v1'
    for key,name in [('compiler_owner','compile-endpoint-normal/owner.json'),('output_seal','compile-endpoint-normal/output-seal.json'),('compiler_products','compiler-products.json')]:joined(R/P['compiler_out']/name,acceptance[key])
    assert not acceptance['structural_graph_accepted'] and not acceptance['runtime_or_numerical_admission']
    terminal=json.loads((old/'independent-terminal-review.json').read_text());assert terminal['passed'] and terminal['source_revision']==P['operator_revision']
    for name,pin in terminal['critical_receipts'].items():joined(R/P['compiler_out']/name,pin)
    source_review=json.loads((old/'independent-source-review.json').read_text());assert source_review['passed']
    for name,review_name in [('root-source-admission.json','independent-source-review.json'),('root-terminal-admission.json','independent-terminal-review.json')]:
        admission=json.loads((old/name).read_text());assert admission['passed'] and admission['source_revision']==P['operator_revision']
        joined(R/P['compiler_out']/'source-freeze.json',admission['source_freeze'])
        joined(R/P['compiler_out']/review_name,admission['independent_review'])
    assert json.loads((old/'root-terminal-admission.json').read_text())['compiler_producer_accepted']
    origins=json.loads((old/'origins.json').read_text());stop=origins['previous_header_and_frozen_builder']
    failed_path=Path(stop['archive'])/'files/out/portable-divide-fp32-estimate-screen/compile-endpoint-normal/owner.json'
    failed=json.loads(failed_path.read_text());assert not failed['passed'] and failed['child_exit']==1 and failed['status']=='terminal'
    assert failed['source_unchanged'] and failed['whole_inputs_unchanged'] and failed['remaining_owned_processes']==[] and failed['cleanup_errors']==[]
    assert terminal['lineage']['original_failed_compiler_owner']['receipt']==record(failed_path)
    assert result['previous_actual_cost_premise']['decision']=='COST_PREMISE_SURVIVES_UNPROVEN'
    for key in ('commit','cost_result','manifest','root_admission','decision'):
        assert terminal['lineage']['previous_distinct_UDiv_correction_JIT'][key]==result['previous_actual_cost_premise'][key]
    return {'compiler_commit':P['compiler_commit'],'compiler_manifest':record(manifest),'compiler_acceptance':record(old/'compiler-producer-acceptance.json'),
            'compiler_owner':record(folder/'owner.json'),'compiler_seal':record(folder/'output-seal.json'),'compiler_products':record(old/'compiler-products.json'),
            'compiler_terminal_review':record(old/'independent-terminal-review.json'),'compiler102_static_and636_consumed_inputs_verified':True,
            'compiler_independent_source_review':record(old/'independent-source-review.json'),
            'compiler_root_source_admission':record(old/'root-source-admission.json'),'compiler_root_terminal_admission':record(old/'root-terminal-admission.json'),
            'explicit_retired_root_relocation_preserves_historical_resolved_identity':True,
            'stock_module':result['modules']['stock']['module'],'stock_assembly':result['modules']['stock']['assembly'],
            'raw_assembly':record(products()/P['product_names'][1]),'raw_module':record(products()/P['product_names'][0]),
            'candidate_products':{name:record(products()/name) for name in P['product_names']},
            'completed_compiler_policy_receipt_reused_without_rerun':True,'stock_source_tool_receipt':result['stock_reuse'],
            'original_failed_owner':record(failed_path),'original_failed_passed':False,'original_failed_child_exit':1,
            'previous_UDiv_correction_cost_premise_survived_unproven':True,
            'fresh_saved_birth_absence':{'accepted_compiler':fresh_births(owner),'original_failed_compiler':fresh_births(failed)}}
