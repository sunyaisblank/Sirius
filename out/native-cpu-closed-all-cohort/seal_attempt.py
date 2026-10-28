"""Seal one original AX frame outcome after fresh all8/22/Decimal; run no gate."""
from collections import Counter
from pathlib import Path
import hashlib
import importlib.util
import json
import subprocess
import zipfile

R=Path.cwd();Q=Path(__file__).resolve().parent;F=Q/'cold-frame'
C8=R/'attestations/qualification/363aa0c2/native-cpu-closed-all-original-c8/MANIFEST.json'
AUTHORITY=R/'out/seal-native-cpu-shared-value-attempt.py'
LIMITS=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152)


def pin(path):
    with path.open('rb') as stream:return dict(bytes=path.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())
def doc(path):return json.loads(path.read_text())
def relative(path):return str(path.relative_to(R))
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
    capsule=doc(C8);c8=doc(Q/'control-8-acceptance.json')
    assert c8['passed'] and c8['all32_readbacks_match_preserved_originals']
    assert pin(Q/'control-8-acceptance.json')==capsule['original_C8_acceptance']
    assert c8['observed_Transport12_implementation']=='F scalar one-row tasks'
    assert c8['optional_math']['selected_simd_stages']==c8['optional_math']['observed_optional_math_stages']==['dopri_phase','endpoint']
    assert c8['optional_math']['all_selected_optional_math_observed'] and c8['optional_math']['genuine_private_FP_scope_observed']
    benefit_path=Q/'benefit/readers/benefit.json';decimal_path=Q/'benefit/readers/decimal-agreement.json'
    controls_path=Q/'controls/readers/controls-1-8-acceptance.json'
    benefit,decimal,controls=doc(benefit_path),doc(decimal_path),doc(controls_path)
    assert benefit['audit_passed'] and benefit['all22_numeric_pass'] and benefit['checks_total']==benefit['checks_passed']==22
    assert benefit['fixed_order']==['a1','b1','b2','a2'] and benefit['internal_repeats_per_slot']==3 and benefit['extra_retries_or_warmups']==0
    assert decimal['agrees_on_all22'] and decimal['all22_pass'] and len(decimal['checks'])==22
    assert controls['all8_original_controls_pass'] and controls['all22_original_benefit_checks_pass']
    assert controls['control8_acceptance']==pin(Q/'control-8-acceptance.json')
    assert controls['benefit']==pin(benefit_path) and controls['decimal']==pin(decimal_path)
    assert controls['reader']==pin(Q/'controls/readers/read_controls.py')
    assert controls['optional_C8_runtime_coverage']==c8['optional_math']
    assert set(controls['remaining7'])=={'control-'+str(i) for i in range(1,8)}
    benefit_actions=doc(Q/'benefit/runtime/actions.json')
    benefit_capsule=R/'attestations/qualification'/benefit_actions['a1']['revision'][:8]/'native-cpu-closed-all-original-benefit/MANIFEST.json'
    sealed_benefit=doc(benefit_capsule)
    assert sealed_benefit['benefit_gate_passed'] and sealed_benefit['execution_revision']==benefit_actions['a1']['revision']
    assert sealed_benefit['benefit_receipt']==row(benefit_path) and sealed_benefit['decimal_receipt']==row(decimal_path)
    assert sealed_benefit['benefit']==benefit and sealed_benefit['decimal']==decimal
    assert sealed_benefit['original_C8_capsule']==row(C8)
    assert sealed_benefit['executed_original_slots']==['a1','b1','b2','a2'] and not sealed_benefit['unexecuted_original_slots']
    assert all(record['passed'] for record in sealed_benefit['owners'].values())
    provider=doc(Q/'benefit/providers.json')['providers']['candidate']
    assert provider['library']==str(R/capsule['library']['path'])
    assert provider['library_identity']=={k:capsule['library'][k] for k in ('bytes','sha256')}
    for key in ('library','ICD'):assert pin(Path(provider[key]))==provider[key+'_identity']
    # Check the already-prepared namespace correction as source only. Never
    # invoke check_acceptance(), a producer, or the terminal reader here.
    spec=importlib.util.spec_from_file_location('AX_seal_frame_source_authority',F/'frame_authority.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    reader_source=module.load('read_result.py')
    assert "current = root / 'out/native-cpu-closed-all-cohort'" in reader_source
    assert "current = root / 'out/native-cpu-dispatch-token-cohort'" not in reader_source
    owners=[(R/'out/native-cpu-closed-all-component/build',R/'out/native-cpu-closed-all-component/actions.json','build'),
        (R/'out/native-cpu-closed-all-provider/provider/configure',R/'out/native-cpu-closed-all-provider/provider/actions.json','configure'),
        (R/'out/native-cpu-closed-all-provider/provider/build',R/'out/native-cpu-closed-all-provider/provider/actions.json','build'),
        (Q/'control-8',Q/'actions.json','control-8')]
    owners += [(Q/'benefit/runtime'/name,Q/'benefit/runtime/actions.json',name) for name in ('a1','b1','b2','a2')]
    owners += [(Q/'controls/runtime'/('control-'+str(i)),Q/'controls/runtime/actions.json','control-'+str(i)) for i in range(1,8)]
    owners += [(F/'frame',F/'actions.json','frame')]
    history={}
    for folder,actions,label in owners:
        path=folder/'owner.json';owner=doc(path);action=doc(actions)[label]
        assert owner['status']=='terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert owner['source']==dict(revision=action['revision'],status='')
        assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
        assert owner['observed_births'] and all(absent(birth) for birth in owner['observed_births'])
        assert path.with_name('inputs-before.json').read_bytes()==path.with_name('inputs-after.json').read_bytes()
        if folder!=F/'frame':assert owner['passed'] and owner['child_exit']==0 and owner['stop_reason'] is None
        history[relative(folder)]=dict(owner=row(path),passed=owner['passed'],child_exit=owner['child_exit'],
            source=owner['source'],argv=owner['argv'],wall_seconds=owner['wall_seconds'],peak_sampled_rss_kib=owner['peak_sampled_rss_kib'],
            stop_reason=owner['stop_reason'],before=row(folder/'inputs-before.json'),after=row(folder/'inputs-after.json'),
            source_unchanged=True,whole_inputs_unchanged=True,all_owned_births_reaped_now=True)
    assert len(history)==16
    for label,record in controls['remaining7'].items():
        folder=Q/'controls/runtime'/label
        assert record['passed'] and record['owner']==pin(folder/'owner.json') and record['XML']==pin(folder/'gtest.xml')
        assert doc(folder/'owner.json')['limits']==doc(Q/'controls/runtime/actions.json')[label]['limits']==LIMITS
    for label in ('a1','b1','b2','a2'):
        assert sealed_benefit['owners'][label]['owner']==row(Q/'benefit/runtime'/label/'owner.json')
        assert benefit['records'][label]['readbacks']==c8['record']['readbacks']
    owner=doc(F/'frame/owner.json');action=doc(F/'actions.json')['frame'];result=doc(F/'result.json')
    assert owner['limits']==action['limits']==LIMITS
    host_path=R/'out/native-cpu-owned-lower-increments/host-build-result.json';host=doc(host_path)
    executable=R/host['frame_executable'];assert host['passed'] and pin(executable)==host['frame_executable_identity']
    assert action['argv']==['/usr/bin/cmake','-DSIRIUS_RENDER_TESTS='+str(executable),'-P',str(R/'tests/operational/vulkan_runtime_profile.cmake')]
    assert action['runtime_provider']==provider['library'] and action['runtime_icd']==provider['ICD']
    assert owner['declared_output_environment']==dict(GTEST_OUTPUT='xml:'+str(F/'frame/gtest.xml'),VK_DRIVER_FILES=provider['ICD'],SIRIUS_VULKAN_DEVICE='0',MESA_SHADER_CACHE_DISABLE='true')
    assert owner['provider_mappings'] and all(item['owned_session']==owner['child_pid'] and item['maps'] and all(provider['library'] in line for line in item['maps']) for item in owner['provider_mappings'])
    assert result['attempted'] and result['accepted']==owner['passed'] and result['owner']==pin(F/'frame/owner.json')
    assert result['controls']==pin(controls_path) and result['benefit']==pin(benefit_path) and result['decimal']==pin(decimal_path)
    assert result['build']==pin(host_path) and result['frame_inputs']==pin(F/'frame-inputs.json')
    assert result['stop_reason']==owner['stop_reason'] and result['child_exit']==owner['child_exit']
    assert result['wall_seconds']==owner['wall_seconds'] and result['peak_sampled_rss_kib']==owner['peak_sampled_rss_kib']
    assert set(result['module_joins'])=={'Transport','Endpoint','DopriPhase'}
    assert all(row['complete_occurrences']==row['executable_complete_occurrences']==row['frame_executable_complete_occurrences']==1 for row in result['module_joins'].values())
    if result['accepted']:
        assert owner['child_exit']==0 and owner['stop_reason'] is None and owner['wall_seconds']<=300
        assert owner['peak_sampled_rss_kib']<=4194304 and owner['minimum_sampled_available_kib']>=2097152
        assert result['XML']==pin(F/'frame/gtest.xml')
    else:
        assert owner['stop_reason'] is not None or owner['child_exit']!=0
        assert result['XML']==(pin(F/'frame/gtest.xml') if (F/'frame/gtest.xml').exists() else None)
    cleanup=doc(F/'disposable-cleanup.json')
    assert cleanup['status']=='complete' and cleanup['canonical_payload_preserved']
    assert cleanup['resource_names_retired']==32 and cleanup['resource_payload_bytes_reclaimed']==0
    # The terminal reader checked live aliases before their separately recorded
    # retirement. Historical owner snapshots stay authoritative; do not demand
    # retired aliases as live inputs or restore resources during sealing.
    revision=action['revision'];assert len(revision)==40
    destination=R/'attestations/qualification'/revision[:8]/'native-cpu-closed-all-original-frame-attempt'
    assert not destination.exists()
    sources=[Path(__file__).resolve(),Q/'seal_benefit.py',AUTHORITY,C8,benefit_capsule,host_path,
        controls_path,benefit_path,decimal_path,Q/'control-8-acceptance.json',Q/'benefit/providers.json']
    sources += [path for parent in (Q/'controls',F) for path in parent.rglob('*') if path.is_file() and '__pycache__' not in path.parts]
    entries={relative(path):pin(path) for path in sorted(set(sources))}
    destination.mkdir(parents=True);archive_path=destination/'source-and-evidence.zip'
    with zipfile.ZipFile(archive_path,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for name in entries:archive.write(R/name,name)
    with zipfile.ZipFile(archive_path) as archive:
        assert len(archive.namelist())==len(entries) and set(archive.namelist())==set(entries)
        for name,expected in entries.items():
            with archive.open(name) as stream:assert hashlib.file_digest(stream,'sha256').hexdigest()==expected['sha256']
            assert archive.getinfo(name).file_size==expected['bytes'] and pin(R/name)==expected
    assert archive_path.stat().st_size<100*1024*1024
    manifest=dict(status='AX_ORIGINAL_FRAME_PASS_BROADER_QUALIFICATION_OPEN' if result['accepted'] else 'AX_ORIGINAL_FRAME_FAILED_STOP',
        custody_revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),execution_revision=revision,
        candidate_revision=capsule['candidate_revision'],owners=history,
        phase_passes=dict(Counter('passed' if row['passed'] else 'failed' for row in history.values())),
        original_C8_capsule=row(C8),original_benefit_capsule=row(benefit_capsule),
        all8_receipt=row(controls_path),benefit_receipt=row(benefit_path),decimal_receipt=row(decimal_path),
        C8_observed_scalar_stages=['transport'],C8_observed_closed_SIMD_stages=['endpoint','dopri_phase'],
        library=capsule['library'],provider=provider,frame_executable=dict(path=host['frame_executable'],**host['frame_executable_identity']),
        final_frame_result=result,final_frame_owner=owner,disposable_cleanup=cleanup,
        corrected_reader_namespace_from_outset=True,reader_source=pin(F/'read_result.py'),frame_authority=pin(F/'frame_authority.py'),
        inherited_sealer_source=row(AUTHORITY),
        bundle=dict(path=relative(archive_path),**pin(archive_path),members=entries),
        source_authority='Exact owner source revisions, argv/limits, before/after frozen snapshots and current immutable capsule/provider/host identities remain authoritative. Compiler/native/library products are recovered from the unchanged C8/context capsules; the sole ABBA is recovered from its fresh benefit capsule. No old receipt, failed outcome, frame input or reader is rewritten.',
        scope='One materially changed AX closed3/local9 implementation, same original C8/all32 and truthful first T12scalar/E-DP SIMD coverage, sole ABBA22 plus independent Decimal, seven original controls/final all8, then one original128x128/1SPP Kerr frame under300s4GiB2GiB with unchanged cache/governors/scenes/layouts/allocations/counts/FP/source tuples and H executable/resources. The original test validates/removes the PPM. Only current actual frame result determines this narrow verdict; a failed frame stays failed and stops downstream acceptance. No retry/extra study/timing transfer/full scientific/physical-hardware/release/main or ultimate-goal acceptance is synthesized.')
    write(destination/'MANIFEST.json',manifest)
    print(json.dumps(dict(sealed=relative(destination),status=manifest['status'],owners=16,
        phase_passes=manifest['phase_passes'],zip_bytes=archive_path.stat().st_size)))


if __name__=='__main__':seal()
