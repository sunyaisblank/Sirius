"""Pinned successor of the retained-artifact audit; only BB comment association differs."""
from pathlib import Path
import hashlib

R=Path.cwd();W=Path(__file__).resolve().parent
OLD=W/'endpoint_residual.py'
OLD_SHA='7890dac08519c8cdfc63ba9135729c51dcb2901f9bf259da617bb659cb6860ee'


def replace(text,before,after):
    assert text.count(before)==1,before
    return text.replace(before,after,1)


HISTORY=r'''
def canonical_history():
    old_cpp=V/'audit_closed_endpoint_residual.cpp'
    checked(old_cpp,dict(bytes=16763,sha256='941c668609446d963e8e60381c4b0fa1343fec75d4a8dc9a460ff9598e7d02f4'))
    checked(W/'endpoint_residual.py',dict(bytes=16070,sha256='7890dac08519c8cdfc63ba9135729c51dcb2901f9bf259da617bb659cb6860ee'))
    terminal(W/'cleanup-audit-plugin-endpoint/owner.json')
    checked(W/'cleanup-audit-endpoint/owner.json',dict(bytes=1711,sha256='e6d5a88b46df5d6c81f7501a8fceac833d04230bc3bed54cf6ba4262d133711f'))
    terminal(W/'cleanup-audit-endpoint/owner.json',False)
    audit_before=doc(W/'cleanup-audit-endpoint/inputs-before.json')
    for old_product in (B/'endpoint.residual-audit.so',W/'endpoint-residual-plugin-result.json',W/'endpoint-residual-plan.json'):
        checked(old_product,audit_before[relative(old_product)])
    checked(B/'endpoint.residual-audit.so',doc(W/'endpoint-residual-plugin-result.json')['output'])
    capture=W/'cleanup-audit-error-capture-endpoint/owner.json';captured=doc(capture)
    assert captured['status']=='terminal' and captured['passed'] is True and captured['child_exit']==1
    assert captured['source_unchanged'] and captured['whole_inputs_unchanged']
    assert not captured['remaining_owned_processes'] and not captured['cleanup_errors']
    assert captured['observed_births'] and all(x['absent'] for x in captured['observed_births'])
    assert capture.with_name('inputs-before.json').read_bytes()==capture.with_name('inputs-after.json').read_bytes()
    original=doc(W/'endpoint-residual-plan.json')
    error_plan=doc(W/'endpoint-residual-error-capture-plan.json')
    assert captured['argv']==error_plan['argv']==original['commands']['audit'] and error_plan['expected_error_exit']==1
    assert 'full canonical module differs from independent retained reference' in capture.with_name('stderr.log').read_text()
    diagnostic=W/'endpoint-residual-serialization-diagnostic'
    terminal(diagnostic/'owner.json')
    checked(diagnostic/'comment-order-proof.json',dict(bytes=723,sha256='f5426c31acec677a5b976c947ce4cc10d6f1b14e1320e4e899d1880594ce8b48'))
    checked(diagnostic/'serialization-comparison.json',dict(bytes=7300,sha256='197a2e5720ad72c5ae3c3a327c136da26057263fcb0c8ef000811f8f50ab950b'))
    proof=doc(diagnostic/'comment-order-proof.json');comparison=doc(diagnostic/'serialization-comparison.json')
    for key in ('passed','only_predecessor_comment_order_differs','all_noncomment_bytes_identical',
                'predecessor_multisets_identical','ignored_module_identifier_only'):assert proof[key] is True,key
    assert proof['reordered_comment_lines']==comparison['differing_lines']==3753
    assert comparison['optimizations_run'] is False and comparison['accepted_artifact'] is False
    recovery=doc(diagnostic/'serialized-text-retirement.json')
    assert recovery['only_generated_comments_reordered'] is True and recovery['comparison_proof']=='comment-order-proof.json'
    assert len(recovery['files'])==2
    # Recovery hashes/proof are live authority. The two diagnostic serializer
    # products are deliberately neither opened nor required nor frozen here.
    for item,filename,source,identity in zip(recovery['files'],('from-bc.ll','from-text.ll'),(BC,LL),(BC_PIN,LL_PIN)):
        assert item['path']==relative(diagnostic/filename) and item['bytes']==119473814
        assert item['recover_from']==dict(path=relative(source),**identity)
    assert recovery['files'][0]['sha256']==LL_PIN['sha256']
    assert recovery['files'][1]['sha256']=='7c8fcac1ccdfdbec214cdc7551c5825cf6e090df56805cff8c56f5a101b72130'
    paths=[old_cpp,W/'endpoint_residual.py',W/'endpoint-residual-plan.json',
        W/'endpoint-residual-prior-actions.json',W/'endpoint-residual-plugin-result.json',
        W/'endpoint-residual-error-capture-plan.json',W/'endpoint-residual-error-capture-prior-actions.json',
        W/'endpoint-residual-serialization-prior-actions.json',B/'endpoint.residual-audit.so',
        diagnostic/'comment-order-proof.json',diagnostic/'serialization-comparison.json',diagnostic/'serialized-text-retirement.json']
    for directory in ('cleanup-audit-plugin-endpoint','cleanup-audit-endpoint',
                      'cleanup-audit-error-capture-endpoint','endpoint-residual-serialization-diagnostic'):
        paths += [W/directory/name for name in ('owner.json','inputs-before.json','inputs-after.json','stdout.log','stderr.log','samples.jsonl')]
    return {relative(p):pin(p) for p in paths}
'''


def prepared_source():
    assert hashlib.sha256(OLD.read_bytes()).hexdigest()==OLD_SHA
    text=OLD.read_text()
    for before,after in (
        ('audit_closed_endpoint_residual.cpp','audit_closed_endpoint_residual_canonical.cpp'),
        ('cleanup-audit-plugin-endpoint','cleanup-canonical-audit-plugin-endpoint'),
        ('cleanup-audit-endpoint','cleanup-canonical-audit-endpoint'),
        ('endpoint-residual-plan.json','endpoint-canonical-residual-plan.json'),
        ('endpoint-residual-proof.json','endpoint-canonical-residual-proof.json'),
        ('endpoint.residual-audit.so','endpoint.residual-canonical-audit.so'),
        ('endpoint-residual-plugin-result.json','endpoint-canonical-residual-plugin-result.json'),
        ('endpoint-residual-prior-actions.json','endpoint-canonical-residual-prior-actions.json'),
        ('sirius-audit-closed-endpoint-residual','sirius-audit-closed-endpoint-residual-canonical')):
        assert before in text,before;text=text.replace(before,after)
    text=replace(text,"    return inputs\n\n\ndef action(",
        "    inputs.update(canonical_history())\n    return inputs\n\n\ndef action(")
    text=replace(text,"    return plan\n\n\ndef freeze_history(",
        "    canonical_history()\n    return plan\n\n\ndef freeze_history(")
    text=replace(text,"    assert proof['numerical_qualification_transferred'] is False",
        "    assert proof['text_reference_predecessor_comments_normalized'] is True\n"+
        "    assert proof['text_program_fingerprint'] and proof['canonical_module_fingerprint']\n"+
        "    assert proof['numerical_qualification_transferred'] is False")
    before="""    original_bytes=(W/'actions.json').read_bytes()
    original_snapshot=doc(W/'cleanup-endpoint/inputs-before.json')[relative(W/'actions.json')]
    assert dict(bytes=len(original_bytes),sha256=hashlib.sha256(original_bytes).hexdigest())=={k:original_snapshot[k] for k in ('bytes','sha256')}"""
    after="""    original_snapshot=doc(W/'cleanup-endpoint/inputs-before.json')[relative(W/'actions.json')]
    checked(W/'endpoint-residual-prior-actions.json',original_snapshot)
    original_actions_subset=doc(W/'endpoint-residual-prior-actions.json')
    current_actions=doc(W/'actions.json')
    assert all(current_actions.get(k)==v for k,v in original_actions_subset.items())
    historical_additions={'cleanup-audit-plugin-endpoint','cleanup-audit-endpoint',
        'cleanup-audit-error-capture-endpoint','endpoint-residual-serialization-diagnostic'}
    assert set(current_actions)==set(original_actions_subset)|historical_additions"""
    text=replace(text,before,after)
    before="""        completed=subprocess.run(plan['commands']['audit'],check=True,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if completed.stderr:print(completed.stderr,file=sys.stderr,end='')
        proof=json.loads(completed.stdout);proof_checked(proof)"""
    after="""        completed=subprocess.run(plan['commands']['audit'],check=False,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if completed.stderr:print(completed.stderr,file=sys.stderr,end='')
        if completed.returncode:
            if completed.stdout:print(completed.stdout,end='')
            raise SystemExit(completed.returncode)
        proof=json.loads(completed.stdout);proof_checked(proof)"""
    text=replace(text,before,after)
    text=replace(text,'            original_failed_owner=FAILED_PIN,original_context_plan=CONTEXT_PLAN_PIN,',
        "            original_failed_owner=FAILED_PIN,original_context_plan=CONTEXT_PLAN_PIN,\n"+
        "            original_failed_audit_owner=row(W/'cleanup-audit-endpoint/owner.json'),\n"+
        "            predecessor_comment_diagnostic=row(W/'endpoint-residual-serialization-diagnostic/comment-order-proof.json'),")
    text=replace(text,'phase,=sys.argv[1:]\n',HISTORY+'\nphase,=sys.argv[1:]\n')
    return text


if __name__=='__main__':
    exec(compile(prepared_source(),str(Path(__file__)),'exec'))
