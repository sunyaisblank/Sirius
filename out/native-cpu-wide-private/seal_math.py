"""Seal only three actual wide-wrapper compiler gates; execute no compiler/workload."""
from pathlib import Path
import hashlib
import json
import subprocess
import zipfile

R=Path.cwd();W=Path(__file__).resolve().parent;RUN=W/'math'
REVISION='d81ace28cdb58445898649a870f014bfdfe77be0'
ORDER=('math-source','math-verify','math-native')
RECIPES=('dense_add','dense_multiply','dense_divide','dense_sqrt','tagged_add','tagged_subtract','tagged_multiply','tagged_divide','tagged_sqrt')
OLD=R/'attestations/qualification/b2cbd4db/native-cpu-private-context-interface-gates/MANIFEST.json'


def pin(path):
    with path.open('rb') as stream:return dict(bytes=path.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())
def identity(data):return dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
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
    assert W==R/'out/native-cpu-wide-private'
    custody_revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    custody_status=subprocess.check_output(['git','status','--porcelain'],text=True)
    assert pin(OLD)==dict(bytes=26762,sha256='47793db110e61869ed36dce41b0386e747068c4f8ca894fead65d7ffffefdf69')
    old=doc(OLD);plan=doc(W/'math-plan.json');actions=doc(RUN/'actions.json')
    assert list(plan['commands'])==list(ORDER) and set(actions)==set(ORDER)
    results={name:doc(W/(name+'-result.json')) for name in ORDER}
    contract=doc(W/'math-source-contract.json')
    assert len(contract['recipes'])==len(contract['exact_old_core_headers'])==9
    assert contract['output']==results['math-source']['output']
    history={};snapshots={};sources={};source_custody={};evidence=[];products=[]
    for index,name in enumerate(ORDER):
        folder=RUN/name;owner=doc(folder/'owner.json');action=actions[name]
        assert owner['status']=='terminal' and owner['passed'] and owner['child_exit']==0
        assert owner['source']==dict(revision=REVISION,status='') and action['revision']==REVISION
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert owner['limits']==action['limits']==dict(seconds=300,rss_kib=4194304,reserve_kib=2097152)
        assert owner['argv']==['/usr/bin/python3','-B',str(W/'build_math.py'),name]
        assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
        assert owner['observed_births'] and all(absent(birth) for birth in owner['observed_births'])
        assert (folder/'inputs-before.json').read_bytes()==(folder/'inputs-after.json').read_bytes()
        before=doc(folder/'inputs-before.json');snapshots[name]=before
        assert {k:before[relative(W/'math-plan.json')][k] for k in ('bytes','sha256')}==pin(W/'math-plan.json')
        assert {k:before[relative(RUN/'actions.json')][k] for k in ('bytes','sha256')}==pin(RUN/'actions.json')
        result=results[name];assert result['passed'] and result['workload'] is False
        output=R/result['output']['path'];assert pin(output)=={k:result['output'][k] for k in ('bytes','sha256')}
        assert output.name=={'math-source':'wide-math.ll','math-verify':'wide-math.bc','math-native':'wide-math.o'}[name]
        products.append(output)
        incoming=None
        if index:
            previous=results[ORDER[index-1]]['output'];entry=before[previous['path']]
            assert {k:entry[k] for k in ('bytes','sha256')}=={k:previous[k] for k in ('bytes','sha256')}
            assert Path(entry['resolved_path']).resolve()==(R/previous['path']).resolve()
            incoming=dict(previous_result=row(W/(ORDER[index-1]+'-result.json')),
                frozen_input=dict(path=previous['path'],**{k:entry[k] for k in ('bytes','sha256')}),exact_join=True)
        history[name]=dict(owner=row(folder/'owner.json'),source=owner['source'],passed=True,
            wall_seconds=owner['wall_seconds'],peak_sampled_rss_kib=owner['peak_sampled_rss_kib'],
            before=row(folder/'inputs-before.json'),after=row(folder/'inputs-after.json'),prior_output_join=incoming,
            whole_inputs_unchanged=True,source_unchanged=True,all_owned_births_reaped_now=True)
        evidence += [path for path in folder.iterdir() if path.is_file()]
        evidence.append(W/(name+'-result.json'))
    old_products={}
    for name,expected in contract['input_math'].items():
        path=R/'bin/linux-gcc/native-cpu-local-effects'/name;assert pin(path)==expected
        key=relative(path);assert plan['inputs'][key]==expected
        archived=[dict(bundle=label,archive={k:bundle[k] for k in ('path','bytes','sha256')},member=key,identity=expected)
            for label,bundle in old['bundles'].items() if bundle.get('members',{}).get(key)==expected]
        assert len(archived)==1
        for stage in ORDER:assert {k:snapshots[stage][key][k] for k in ('bytes','sha256')}==expected
        old_products[key]=dict(identity=expected,recovery=archived[0],unchanged=True)
    for stage in ('math-verify','math-native'):
        key=relative(W/'math-source-contract.json')
        assert {k:snapshots[stage][key][k] for k in ('bytes','sha256')}==pin(W/'math-source-contract.json')
    native=results['math-native'];assert native['output']==dict(path='bin/linux-gcc/native-cpu-wide-private/wide-math.o',bytes=7400,
        sha256='919bd088a551000519a1f358d4886db7a08cea801e0319fb47e216d0203b8b3c')
    wanted={'sirius_wide_value_'+recipe+'16' for recipe in RECIPES}
    assert set(native['native'])==wanted
    references=0
    for recipe in RECIPES:
        name='sirius_wide_value_'+recipe+'16';check=native['native'][name]
        assert check['hidden_regcall_value_entry']==name and check['instructions']>0 and not check['emitted_fused_instructions']
        assert check['core_relocations']=={'sirius_local_value_'+recipe+'8':2};references+=2
        table=check['unwind_table'];assert table['full_function_FDE'] and table['covered_call_sites']==2
        assert table['relative_range']==[0,table['function_bytes']] and table['relative_CFI_rows'] and table['checked_return_offsets']
    assert references==18
    # Source revisions can advance independently after these completed owners.
    # Archive precisely their frozen source image, not an unexecuted edit.
    for name,expected in plan['inputs'].items():
        path=R/name
        if path.suffix not in ('.py','.json') or not path.is_relative_to(R):continue
        for stage in ORDER:assert {k:snapshots[stage][name][k] for k in ('bytes','sha256')}==expected
        current=path.read_bytes() if path.is_file() else None
        if current is not None and identity(current)==expected:
            data=current;method='current_bytes_match_frozen_source'
        else:
            data=subprocess.check_output(['git','show',REVISION+':'+name]);method='exact_execution_revision_Git_blob'
        assert identity(data)==expected,(name,'Original frozen source must be recoverable exactly')
        sources[name]=data
        source_custody[name]=dict(executed_identity=expected,recovered_from=method,execution_revision=REVISION,
            current_identity=identity(current) if current is not None else None)
    evidence += [W/'math-plan.json',RUN/'actions.json',W/'math-source-contract.json',OLD,Path(__file__).resolve()]
    source_entries={relative(path):pin(path) for path in sorted(set(evidence))}
    assert not set(sources)&set(source_entries)
    destination=R/'attestations/qualification'/REVISION[:8]/'native-cpu-wide-private-math-wrapper-gates'
    assert not destination.exists();destination.mkdir(parents=True)
    bundles={}
    for label,paths,originals in (('source-and-evidence',source_entries,sources),
            ('products',{relative(path):pin(path) for path in products},{})):
        entries=dict(paths);entries.update({name:identity(data) for name,data in originals.items()})
        target=destination/(label+'.zip')
        with zipfile.ZipFile(target,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
            for name in paths:archive.write(R/name,name)
            for name,data in originals.items():archive.writestr(name,data)
        with zipfile.ZipFile(target) as archive:
            assert len(archive.namelist())==len(entries) and set(archive.namelist())==set(entries)
            for name,expected in entries.items():
                with archive.open(name) as stream:assert hashlib.file_digest(stream,'sha256').hexdigest()==expected['sha256']
                assert archive.getinfo(name).file_size==expected['bytes']
                if name not in originals:assert pin(R/name)==expected
        assert target.stat().st_size<100*1024*1024
        bundles[label]=dict(path=relative(target),**pin(target),members=entries)
    manifest=dict(status='WIDE_VALUE_WRAPPER_SOURCE_VERIFY_NATIVE_OBJECT_PASS_RUNTIME_UNQUALIFIED',
        execution_revision=REVISION,custody_revision=custody_revision,custody_git_status=custody_status,
        owners=history,phase_passes=dict(passed=3,failed=0),commands=plan['commands'],results=results,
        native_object=native['output'],wrapper_count=9,old_core_relocation_references=18,full_function_FDEs_recorded=9,
        emitted_FMA_instructions=0,source_custody=source_custody,original_math_capsule=row(OLD),
        unchanged_old_math_products=old_products,bundles=bundles,
        numerical_qualification_transferred=False,native_ABI_runtime_qualified=False,provider_adopted=False,
        source_authority='Each later owner frozen-input identity equals the prior result output. Exact executed source bytes must match all three original snapshots; a changed current tracked source is recovered from the d81ace28 Git blob and checked against that frozen identity. No old receipt/source is rewritten.',
        scope='Three completed source-derivation/LLVM-verification/native-wrapper compiler gates only: nine integer-only value16 wrappers, eighteen named relocation references to unchanged strict value8 cores, no emitted FMA and recorded whole-function FDE coverage. No final linked target/complete saved-GP/native ABI or arithmetic execution, raw16 subgroup/context/dispatch, current numerical C8, performance, frame, provider/scientific/release acceptance is claimed or transferred. Original expensive math products are referenced by exact existing recovery-capsule identities, not recopied.')
    write(destination/'MANIFEST.json',manifest)
    print(json.dumps(dict(sealed=relative(destination),owners=3,wrapper_functions=9,old_core_references=18,
        bundles={label:value['bytes'] for label,value in bundles.items()},workload=False)))


if __name__=='__main__':seal()
