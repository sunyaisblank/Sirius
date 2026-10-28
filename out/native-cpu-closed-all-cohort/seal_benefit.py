"""Seal the sole AX ABBA outcome; no reader, workload or broader gate is run."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import subprocess
import zipfile

R=Path.cwd();Q=Path(__file__).resolve().parent;W=Q/'benefit';RUN=W/'runtime';READ=W/'readers'
ORDER=('a1','b1','b2','a2');LIMITS=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152)
C8=R/'attestations/qualification/363aa0c2/native-cpu-closed-all-original-c8/MANIFEST.json'
AUTHORITY=R/'out/seal-native-cpu-shared-value-attempt.py'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())

def doc(path):return json.loads(path.read_text())
def relative(path):return str(path.relative_to(R)) if path.is_relative_to(R) else str(path)
def row(path):return dict(path=relative(path),**pin(path))

def write(path,value):
    with path.open('x') as stream:json.dump(value,stream,indent=2);stream.write('\n')


def absent(birth):
    try:current=(Path('/proc')/str(birth['pid'])/'stat').read_text().rsplit(')',1)[1].split()[19]
    except FileNotFoundError:return True
    return current!=birth['saved_start_ticks']


def seal():
    assert Q==R/'out/native-cpu-closed-all-cohort'
    assert not subprocess.check_output(['git','status','--porcelain']).strip()
    assert pin(AUTHORITY)==dict(bytes=7535,sha256='d4797e1d4242342c2ba98bdf8b9b8003ad416d104360731de89da9841fa3fd9c')
    assert pin(C8)['sha256']=='7603fbd5ed5f98178550ff78ac04b4b64603aecd70d14cd29a42fb1e2eb397f6'
    capsule=doc(C8);accepted=doc(Q/'control-8-acceptance.json')
    assert pin(Q/'control-8-acceptance.json')==capsule['original_C8_acceptance']
    assert accepted['passed'] and accepted['all32_readbacks_match_preserved_originals']
    assert accepted['observed_Transport12_implementation']=='F scalar one-row tasks'
    assert accepted['optional_math']['selected_simd_stages']==accepted['optional_math']['observed_optional_math_stages']==['dopri_phase','endpoint']
    assert accepted['optional_math']['all_selected_optional_math_observed'] and accepted['optional_math']['genuine_private_FP_scope_observed']
    assert capsule['C8_observed_scalar_stages']==['transport'] and capsule['C8_observed_closed_SIMD_stages']==['endpoint','dopri_phase']
    actions=doc(RUN/'actions.json');providers=doc(W/'providers.json')['providers']
    original=doc(R/'out/native-original-input-ownership/runtime/actions.json')
    revision=actions['a1']['revision']
    assert len(revision)==40 and all(c in '0123456789abcdef' for c in revision)
    for label in ORDER:
        action=actions[label]
        for key in ('arm','case','prerequisites','kind','expected_exit','limits'):assert action[key]==original[label][key],(label,key)
        assert action['revision']==revision and action['limits']==LIMITS
        assert action['argv']==['/usr/bin/python3','-B',str(RUN/'run_case.py'),label,'{folder}']
        provider=providers[action['arm']]
        assert action['runtime_provider']==provider['library'] and action['runtime_icd']==provider['ICD']
        assert action['software_provider'] is True
        for key in ('library','ICD'):assert pin(Path(provider[key]))==provider[key+'_identity']
    assert providers['baseline']['library_identity']['sha256']=='cde2d0f8f243448a3f3b28408121412fb03a8bd9c01f5731f9e87a107150e952'
    assert providers['candidate']['library']==str(R/capsule['library']['path'])
    assert providers['candidate']['library_identity']=={k:capsule['library'][k] for k in ('bytes','sha256')}
    assert pin(Path(providers['candidate']['library']))==providers['candidate']['library_identity']
    host=doc(R/'out/native-cpu-owned-lower-increments/host-build-result.json')
    assert host['passed'] and all(action['executable']==str(R/host['executable']) for label,action in actions.items() if label in ORDER)
    assert pin(R/host['executable'])==host['executable_identity']
    sole=doc(READ/'sole-ABBA-inputs.json')
    assert sole['order']==list(ORDER) and sole['extra_warmups_or_resampling']==0
    cache={}
    def current_pin(path):
        key=str(path.resolve())
        if key not in cache:cache[key]=pin(path) if path.is_file() else None
        return cache[key]
    for name,expected in sole['inputs'].items():assert current_pin(R/name)==expected,name
    launched=[label for label in ORDER if (RUN/label).exists()]
    assert launched and launched==list(ORDER[:len(launched)]), 'Only the original ordered executed prefix can be sealed.'
    history={};owners={};custody_failures=[];last=None;births=set()
    for label in launched:
        folder=RUN/label;path=folder/'owner.json';owner=doc(path);action=actions[label]
        assert owner['status']=='terminal'
        assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
        assert owner['observed_births'] and all(absent(birth) for birth in owner['observed_births'])
        before=doc(folder/'inputs-before.json');after=doc(folder/'inputs-after.json')
        equal=(folder/'inputs-before.json').read_bytes()==(folder/'inputs-after.json').read_bytes()
        expected_source=dict(revision=action['revision'],status='')
        source_equal=owner['source']==expected_source
        expected_argv=[word.replace('{folder}',str(folder)) for word in action['argv']]
        argv_equal=owner['argv']==expected_argv
        assert owner['limits']==action['limits']==LIMITS
        changed=[]
        for name,identity in before.items():
            identity_pin={k:identity[k] for k in ('bytes','sha256')}
            if current_pin(R/name)!=identity_pin:changed.append(name)
        identities=all(current_pin(R/name)==identity for name,identity in action['expected_identities'].items())
        birth=(owner['boot_id'],owner['child_pid'],owner['child_start_ticks'])
        assert birth not in births;births.add(birth)
        if last is not None:assert birth[0]==last[0] and int(birth[2])>int(last[2])
        last=birth
        clean=equal and source_equal and argv_equal and not changed and identities and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        if not clean:custody_failures.append(label)
        if owner['passed']:
            assert owner['child_exit']==0 and owner['stop_reason'] is None
            assert owner['wall_seconds']<=300 and owner['peak_sampled_rss_kib']<=4194304
            assert owner['minimum_sampled_available_kib']>=2097152
        history[label]=dict(owner=row(path),passed=owner['passed'],child_exit=owner['child_exit'],source=owner['source'],
            expected_source=expected_source,argv=owner['argv'],expected_argv=expected_argv,
            source_unchanged=owner['source_unchanged'],whole_inputs_unchanged=owner['whole_inputs_unchanged'],
            exact_before_after_equal=equal,current_frozen_input_mismatches=changed,expected_identities_match=identities,
            before=row(folder/'inputs-before.json'),after=row(folder/'inputs-after.json'),
            wall_seconds=owner['wall_seconds'],peak_sampled_rss_kib=owner['peak_sampled_rss_kib'],
            minimum_sampled_available_kib=owner['minimum_sampled_available_kib'],stop_reason=owner['stop_reason'],
            all_owned_births_reaped_now=True)
        owners[label]=owner
    failed=[label for label in launched if not owners[label]['passed']]
    if failed:
        assert failed==[launched[-1]], 'Original runner stops at its first failed owner.'
    if len(launched)<4:assert failed, 'No invented completion of unexecuted original slots.'
    benefit_path=READ/'benefit.json';decimal_path=READ/'decimal-agreement.json'
    benefit=doc(benefit_path) if benefit_path.exists() else None
    decimal=doc(decimal_path) if decimal_path.exists() else None
    if benefit is not None:
        assert len(launched)==4 and not failed and not custody_failures
        assert benefit['audit_passed'] and benefit['checks_total']==22
        assert benefit['fixed_order']==list(ORDER) and benefit['internal_repeats_per_slot']==3 and benefit['extra_retries_or_warmups']==0
        checks=benefit['checks']+benefit['accounting_checks'];assert len(checks)==22
        assert benefit['checks_passed']==sum(check['pass'] for check in checks)
        assert benefit['passed']==benefit['benefit_passed']==all(check['pass'] for check in checks)
        assert set(benefit['records'])==set(benefit['owners'])==set(ORDER)
        assert len(accepted['record']['readbacks'])==32
        for label in ORDER:
            assert benefit['records'][label]['readbacks']==accepted['record']['readbacks']
            assert benefit['owners'][label]['owner']==pin(RUN/label/'owner.json')
    if decimal is not None:
        assert benefit is not None and benefit['audit_passed']
        assert decimal['agrees_on_all22'] and len(decimal['checks'])==22
        assert decimal['all22_pass']==all(check['passed'] for check in decimal['checks'])==benefit['benefit_passed']
        normal={(x['direction'],x.get('family','instance'),x['scope']):x['pass'] for x in benefit['checks']+benefit['accounting_checks']}
        exact={(x['direction'],x['family'],x['scope']):x['passed'] for x in decimal['checks']}
        assert normal==exact and len(exact)==22
    passed=not failed and not custody_failures and len(launched)==4 and benefit is not None and benefit['benefit_passed'] and decimal is not None and decimal['agrees_on_all22'] and decimal['all22_pass']
    status=('AX_SOLE_ORIGINAL_ABBA22_DECIMAL_PASS_BROADER_GATES_OPEN' if passed else
        'AX_SOLE_ORIGINAL_ABBA_OWNER_FAILED_STOP' if failed else
        'AX_SOLE_ORIGINAL_ABBA_CUSTODY_FAILED_STOP' if custody_failures else
        'AX_SOLE_ORIGINAL_ABBA22_NUMERIC_FAILED_STOP' if benefit is not None and not benefit['benefit_passed'] else
        'AX_SOLE_ORIGINAL_ABBA_READERS_UNACCEPTED_STOP')
    destination=R/'attestations/qualification'/revision[:8]/'native-cpu-closed-all-original-benefit'
    assert not destination.exists()
    paths=[Path(__file__).resolve(),AUTHORITY,C8,Q/'control-8-acceptance.json',Q/'control-8/owner.json',
        Q/'control-8/gtest.xml',Q/'control-8/native-observation.json',Q/'control-8/inputs-before.json',Q/'control-8/inputs-after.json',
        R/'out/native-cpu-closed-all-provider/provider/source-contract.json',
        R/'out/native-cpu-closed-all-provider/provider/provider-build-result.json',
        R/'out/native-cpu-closed-all-provider/provider/linked-native-proof.json',
        R/'out/native-cpu-owned-lower-increments/host-build-result.json']
    paths += [path for path in W.rglob('*') if path.is_file() and '__pycache__' not in path.parts]
    entries={relative(path):pin(path) for path in sorted(set(paths))}
    destination.mkdir(parents=True)
    target=destination/'source-and-evidence.zip'
    with zipfile.ZipFile(target,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for name in entries:archive.write(R/name,name)
    with zipfile.ZipFile(target) as archive:
        assert len(archive.namelist())==len(entries) and set(archive.namelist())==set(entries)
        for name,expected in entries.items():
            with archive.open(name) as stream:assert hashlib.file_digest(stream,'sha256').hexdigest()==expected['sha256']
            assert archive.getinfo(name).file_size==expected['bytes'] and pin(R/name)==expected
    assert target.stat().st_size<100*1024*1024
    manifest=dict(status=status,benefit_gate_passed=bool(passed),custody_revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        execution_revision=revision,owners=history,phase_passes=dict(Counter('passed' if x['passed'] else 'failed' for x in history.values())),
        fixed_order=list(ORDER),executed_original_slots=launched,unexecuted_original_slots=list(ORDER[len(launched):]),
        internal_repeats_per_slot=3,extra_trials_warmups_or_resampling=0,original_resource_limits=LIMITS,
        original_C8_capsule=row(C8),original_C8_acceptance=row(Q/'control-8-acceptance.json'),
        C8_observed_scalar_stages=['transport'],C8_observed_closed_SIMD_stages=['endpoint','dopri_phase'],
        C8_closed_Transport_SIMD_observation_claimed=False,providers=providers,
        executable=dict(path=host['executable'],**host['executable_identity']),sole_ABBA_input_ledger=row(READ/'sole-ABBA-inputs.json'),
        benefit=benefit,decimal=decimal,
        benefit_receipt=row(benefit_path) if benefit is not None else None,
        decimal_receipt=row(decimal_path) if decimal is not None else None,
        all32_original_full_readbacks_verified=benefit is not None and benefit['audit_passed'],
        current_frozen_custody_failures=custody_failures,
        bundle=dict(path=relative(target),**pin(target),members=entries),
        inherited_sealer_source=row(AUTHORITY),
        source_authority='Every actual owner source revision/argv, original before/after freeze, current provider/host identities and unchanged original action protocol remain explicit. Failed owners/results and absent downstream receipts remain failed/absent; no reader, workload, retry or acceptance is manufactured by sealing.',
        scope='One sole original AX A1/B1/B2/A2 benefit attempt on the unchanged H executable/stock A and current closed3/local9 B, same22 predicates, full32 buffers, all original lexemes/Decimal gates and resource/cache/governor/cost accounting. Current C8 is linked from its immutable capsule and honestly observed Transport12 scalar plus E/DP SIMD. No C8 timing or prior AU result transfers benefit. Failed numeric/Decimal/owner outcomes stop downstream acceptance. Even all22/Decimal PASS leaves original remaining controls/final all8, frame, scientific/hardware-domain/release/main acceptance and ultimate goal open.')
    write(destination/'MANIFEST.json',manifest)
    print(json.dumps(dict(sealed=relative(destination),status=status,owners=len(history),
        passed=bool(passed),failed_owners=failed,missing_original_slots=manifest['unexecuted_original_slots'],zip_bytes=target.stat().st_size)))


if __name__=='__main__':seal()
