"""Read immutable original stop and separate successful policy-finisher receipts."""
from pathlib import Path
import ast,hashlib,json,os,stat

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

def archive():return Path(P['recovery_archive'])
def oldwork():return archive()/'files'/P['recovery_out']
def products():return archive()/'files'/P['recovery_bin']

def original_path(name):
    path=Path(name);return R/path if not path.is_absolute() else path

def relocated(name):
    path=original_path(name)
    return archive()/'files'/path.relative_to(R) if any(path.is_relative_to(R/P[key]) for key in ('recovery_out','recovery_bin')) else path

def joined(name,pin):
    original=original_path(name);actual=relocated(name);seen=record(actual)
    if actual!=original:
        # Only retired owned regular-file roots relocate. Raw historical resolved identities remain exact.
        assert not actual.is_symlink() and pin['resolved_path']==str(original),name
        seen['resolved_path']=str(original)
    assert {k:seen[k] for k in pin}==pin,name
    return actual

def historical_record(name):
    original=original_path(name);result=record(relocated(name))
    if relocated(name)!=original:result['resolved_path']=str(original)
    return result

def fresh_births(owner):
    current_boot=Path('/proc/sys/kernel/random/boot_id').read_text().strip();rows=[]
    for birth in owner['observed_births']:
        assert birth['absent']
        path=Path('/proc')/str(birth['pid'])/'stat'
        try:
            text=path.read_text();fields=text[text.rfind(')')+2:].split();ticks=fields[19]
        except (FileNotFoundError,ProcessLookupError):ticks=None
        same=current_boot==owner['boot_id'] and ticks==birth['saved_start_ticks']
        assert not same,('saved birth active',birth)
        rows.append({'pid':birth['pid'],'saved_start_ticks':birth['saved_start_ticks'],'current_start_ticks':ticks,'same_boot':current_boot==owner['boot_id'],'saved_birth_absent':True})
    return rows

def completed():
    a=archive();manifest=a/'manifest.json'
    assert content(manifest)==P['recovery_manifest']
    files=json.loads(manifest.read_text())['files'];assert len(files)==38
    for name,pin in files.items():
        path=a/name;assert not path.is_symlink() and content(path)==pin,name
    old=oldwork();params=json.loads((old/'packet.json').read_text())
    assert params['stop_commit']==P['stop_commit'] and params['stop_manifest_sha256']==P['stop_manifest_sha256']
    # Definition-only reuse of the exact accepted finisher's completed-stop verifier.
    # It verifies all old96/630 records, original false/1, 8 source+16 tool stock reuse, and raw outputs.
    tree=ast.parse((old/'recover.py').read_text())
    code=ast.Module(body=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name=='completed'],type_ignores=[])
    scope={'Path':Path,'P':params,'R':R,'record':record,'content':content,'json':json,'re':__import__('re')}
    exec(compile(code,'immutable_completed_stop_definition','exec'),scope)
    stop=scope['completed']()
    frozen=json.loads((old/'source-freeze.json').read_text())['pins'];assert len(frozen)==662
    for name,pin in frozen.items():joined(name,pin)
    context=json.loads((old/'operator-context.json').read_text());assert context['revision']==P['operator_revision']
    for key,name in [('actions','actions.json'),('source_freeze','source-freeze.json'),('preparation','preparation.json')]:joined(R/P['recovery_out']/name,context[key])
    assert (old/'actions.json').read_bytes()==(old/'actions.unpinned.json').read_bytes()
    owner=json.loads((old/'finish-policy/owner.json').read_text());action=json.loads((old/'actions.json').read_text())['finish-policy']
    assert owner['status']=='terminal' and owner['passed'] and owner['child_exit']==0
    assert owner['stop_reason'] is None and owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[]
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['source']=={'revision':P['operator_revision'],'status':''} and owner['limits']==P['limits']
    assert owner['argv']==action['argv'] and len(owner['observed_births'])==1
    assert (old/'finish-policy/inputs-before.json').read_bytes()==(old/'finish-policy/inputs-after.json').read_bytes()
    before=json.loads((old/'finish-policy/inputs-before.json').read_text());assert len(before)==666
    for name,pin in before.items():joined(name,pin)
    seal=json.loads((old/'finish-policy/output-seal.json').read_text());result=json.loads((old/'policy-result.json').read_text())
    for key,name in [('source_freeze','source-freeze.json'),('operator_context','operator-context.json'),('inputs_before','finish-policy/inputs-before.json'),('result','policy-result.json')]:joined(R/P['recovery_out']/name,seal[key])
    expected={str(R/P['recovery_bin']/name):historical_record(R/P['recovery_bin']/name) for name in P['product_names']}
    assert seal['outputs']==result['outputs']==expected and seal['structural_graph_accepted'] is False
    for name,pin in expected.items():joined(name,pin)
    assert result['passed'] and result['completed_stop']==stop
    assert result['one_annotation_and_original_two_RTE_lines_reversal_equal_raw'] and result['unchanged_frozen_private_and_original_normal_policies_passed']
    assert result['complete_type_interface_execution_modes_and_ordered_barriers_equal_stock'] and result['shared_bound_bytes']==12596
    assert result['old_failed_passed'] is False and result['old_failed_child_exit']==1
    assert not result['compiler_or_SSA_or_old_disassembly_repeated'] and not result['structural_graph_accepted'] and not result['runtime_or_numerical_admission'] and not result['benefit_established']
    acceptance=json.loads((old/'policy-acceptance.json').read_text());assert acceptance['passed'] and acceptance['schema']=='sirius-fp32-estimate-policy-recovery-acceptance-v1'
    for key,name in [('reader_owner','finish-policy/owner.json'),('output_seal','finish-policy/output-seal.json'),('policy_result','policy-result.json')]:joined(R/P['recovery_out']/name,acceptance[key])
    assert acceptance['original_failed_owner']==stop['old_failed_owner'] and acceptance['original_failed_passed'] is False and acceptance['original_failed_child_exit']==1
    assert acceptance['remaining_assembly_validation_complete'] and not acceptance['structural_graph_accepted'] and not acceptance['runtime_or_numerical_admission']
    terminal=json.loads((old/'independent-terminal-review.json').read_text());assert terminal['passed'] and terminal['source_revision']==P['operator_revision']
    for name,pin in terminal['critical_pins'].items():joined(R/P['recovery_out']/name,pin)
    failed=json.loads(Path(stop['old_failed_owner']['resolved_path']).read_text())
    return {'recovery_commit':P['recovery_commit'],'recovery_manifest':record(manifest),'policy_acceptance':record(old/'policy-acceptance.json'),
            'policy_owner':record(old/'finish-policy/owner.json'),'policy_seal':record(old/'finish-policy/output-seal.json'),
            'policy_result':record(old/'policy-result.json'),'policy_terminal_review':record(old/'independent-terminal-review.json'),
            'old_failed_owner':stop['old_failed_owner'],'original_failed_passed':False,'original_failed_child_exit':1,
            'old96_static_and630_consumed_inputs_verified':True,'finisher662_static_and666_consumed_inputs_verified':True,
            'explicit_retired_root_relocation_preserves_historical_resolved_identity':True,
            'stock_module':stop['stock_module'],'stock_assembly':stop['stock_assembly'],
            'frozen_builder':stop['frozen_builder'],'raw_assembly':stop['raw_assembly'],'raw_module':stop['raw_module'],
            'candidate_products':{name:record(products()/name) for name in P['product_names']},
            'fresh_saved_birth_absence':{'original_failed':fresh_births(failed),'accepted_finisher':fresh_births(owner)}}
