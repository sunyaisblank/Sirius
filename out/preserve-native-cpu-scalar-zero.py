"""Preserve U's original C3 timeout after the genuine finite/ABBA passes."""
from pathlib import Path
import ast
import hashlib
import json
import stat
import subprocess
import zipfile

R=Path.cwd();W=R/'out/native-cpu-scalar-zero';C=R/'out/native-cpu-scalar-zero-cohort';D=R/'out/native-cpu-scalar-zero-compatibility'
assert not subprocess.check_output(['git','status','--porcelain']).strip()
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
basis=R/'out/preserve-native-cpu-zero-term.py';tree=ast.parse(basis.read_text())
functions=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in ('doc','identity','archive')]
assert {node.name for node in functions}=={'doc','identity','archive'}
exec(compile(ast.Module(body=functions,type_ignores=[]),str(basis),'exec'))
assert doc(C/'control-8-acceptance.json')['passed'] and doc(D/'result.json')['passed']
assert doc(C/'benefit/readers/benefit.json')['all22_numeric_pass']
assert doc(C/'benefit/readers/decimal-agreement.json')['all22_pass']
failed=C/'controls/runtime/control-3/owner.json';timeout=doc(failed)
assert timeout['status']=='terminal' and not timeout['passed'] and timeout['stop_reason']=='diagnostic_timeout'
assert timeout['source']['revision']==revision
assert all(doc(C/'controls/runtime'/('control-'+str(i))/'owner.json')['passed'] for i in (1,2))
assert not any((C/'controls/runtime'/('control-'+str(i))).exists() for i in range(4,8))
assert not (C/'cold-frame/actions.json').exists() and not (C/'cold-frame/resource-restore.json').exists()
owners={}
for root in (W,C,D):
    for path in sorted(root.rglob('owner.json')):
        owner=doc(path)
        assert owner['status']=='terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        assert owner['passed'] or path==failed
        owners[str(path.relative_to(R))]={key:owner.get(key) for key in ('source','status','passed','child_exit','wall_seconds','peak_sampled_rss_kib','stop_reason')}
assert len(owners)==13
outcome={'accepted':False,'actual_new_scalar_objects':3,'new_ISPC_compilations':0,
    'full38_187960_word_compatibility_pass':True,'original_C8_all32_pass':True,
    'sole_ABBA22_and_independent_Decimal_pass':True,'original_controls_1_2_pass':True,
    'original_control3_timeout':timeout['wall_seconds'],
    'unexecuted_after_failure':['control-4','control-5','control-6','control-7','originalframe'],
    'source_revision':revision,
    'reason':'ProjectedEndpointsKeepPhysicalColumnsAndRetainedContinuation hit unchanged300s governor. No rerun or frame; next step is bounded failure-stage isolation, not an assumed arithmetic diagnosis.'}
with (C/'controls/trial-terminal.json').open('x') as stream:
    json.dump(outcome,stream,indent=2);stream.write('\n')
destination=R/'attestations/qualification'/revision[:8]/'native-cpu-scalar-zero-control3-timeout'
assert not destination.exists();destination.mkdir(parents=True)
sources={str(path.relative_to(R)):path for root in (W,C,D) for path in root.rglob('*') if path.is_file()}
for path in (R/'bin/linux-gcc/native-cpu-scalar-zero/typed-batch').iterdir():
    assert path.is_file();sources[str(path.relative_to(R))]=path
for path in (Path(__file__).resolve(),basis):sources[str(path.relative_to(R))]=path
compiled=doc(W/'cpp-compile-result.json');provider=doc(W/'provider/provider-build-result.json');component=doc(D/'build-result.json')
products=[R/row['object']['path'] for row in compiled['objects']]
products += [R/provider[key] for key in ('library','ICD')]
products += [R/component[key]['path'] for key in ('shim_object','executable')]
assert len(products)==len(set(products))==7
bundles=[archive('source-build-and-trial.zip',sources),archive('actual-products.zip',{str(path.relative_to(R)):path for path in products})]
manifest={'schema':'sirius-native-cpu-scalar-zero-terminal-v1','qualifying_source_revision':revision,
    'outcome':outcome,'owners':owners,'archives':bundles,'actual_new_product_members':7,
    'scope':'U applies only the reviewed exact-zero arm and explicit Term inlining to actual D scalar arithmetic, preserving every other header/wrapper byte/all9380calls. Three genuine strict CPP builds, original Q/S scheduling/raw3/admission/adapter reused; exact898function/module joins, fresh full38/all187960words, original C8/all32/51markers/54joins and sole ABBA22+Decimal PASS. Original C1/C2 passed their GTest owned actions; C3 timed out before all8 reader/frame. No independent all-domain math, full scientific/release/main/adoption claim.',
    'additional_qualification_claimed':False,'originalframe_attempted':False,
    'reused_native_objects':doc(W/'provider/source-contract.json')['native_objects'][3:],
    'baseline_driver_runtime_recovery':'attestations/qualification/603a155e/native-cpu-accumulate-inline-frame-timeout/MANIFEST.json',
    'unchanged_S_row_SIMD_recovery':'attestations/qualification/115b86f3/native-cpu-zero-term-frame-timeout/MANIFEST.json',
    'runtime_transition':doc(C/'runtime-transition.json'),
    'additional_profile_executed_in_this_cohort':False}
with (destination/'MANIFEST.json').open('x') as stream:
    json.dump(manifest,stream,indent=2);stream.write('\n')
print(json.dumps({'bundle':str(destination.relative_to(R)),'owners':len(owners),
    'archives':[{'file':row['file'],'bytes':row['identity']['bytes'],'sha256':row['identity']['sha256'],'members':len(row['members'])} for row in bundles]}))
