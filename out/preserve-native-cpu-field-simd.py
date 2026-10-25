"""Preserve the terminal field-SIMD failure; reuse verified ZIP integrity helpers."""
from pathlib import Path
import ast
import hashlib
import json
import stat
import subprocess
import zipfile

R=Path.cwd()
W=R/'out/native-cpu-field-simd'
C=R/'out/native-cpu-field-cohort'
D=R/'out/native-cpu-field-compatibility'
assert not subprocess.check_output(['git','status','--porcelain']).strip()
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
# Load only these already exercised pure archive/identity functions. The old
# preservation workflow and its assertions/actions are never executed.
basis=R/'out/preserve-native-cpu-zero-term.py'
tree=ast.parse(basis.read_text())
functions=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in ('doc','identity','archive')]
assert {node.name for node in functions}=={'doc','identity','archive'}
exec(compile(ast.Module(body=functions,type_ignores=[]),str(basis),'exec'))
assert doc(C/'control-8-acceptance.json')['passed'] and doc(D/'result.json')['passed']
assert doc(C/'control-8-acceptance.json')['field_route_proof']['first_field_layer_returns']==3
failed=C/'benefit/runtime/b1/owner.json';timeout=doc(failed)
assert timeout['status']=='terminal' and not timeout['passed'] and timeout['stop_reason']=='diagnostic_timeout'
assert timeout['source']['revision']==revision
assert doc(C/'benefit/runtime/a1/owner.json')['passed']
assert not any((C/'benefit/runtime'/label).exists() for label in ('b2','a2'))
assert not (C/'controls/runtime/actions.json').exists() and not (C/'cold-frame/actions.json').exists()
assert not (C/'benefit/readers/benefit.json').exists()
owners={}
for root in (W,C,D):
    for path in sorted(root.rglob('owner.json')):
        owner=doc(path)
        assert owner['status']=='terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        if not owner['passed']:
            assert path in (failed,D/'build/owner.json')
        owners[str(path.relative_to(R))]={key:owner.get(key) for key in ('source','status','passed','child_exit','wall_seconds','peak_sampled_rss_kib','stop_reason')}
assert len(owners)==11
partial=doc(D/'build-result.json')
assert len(partial['commands'])==2 and all(row['exit']==0 for row in partial['commands'])
assert doc(D/'build-completion.json')['passed']
outcome={'accepted':False,'original_C8_all32_pass':True,'full38_187960_word_compatibility_pass':True,
    'sole_ABBA_started':['a1','b1'],'baseline_a1_pass':True,'candidate_b1_timeout':timeout['wall_seconds'],
    'all22_and_Decimal_admitted':False,'unexecuted_after_failure':['b2','a2','remaining7controls','originalframe'],
    'reason':'First original candidate benefit slot hit unchanged300s governor; reject field-SIMD T. No changed settings/counts/precision/budgets/cache, no timed rerun or frame attempt.',
    'source_revision':revision}
with (C/'benefit/trial-terminal.json').open('x') as stream:
    json.dump(outcome,stream,indent=2);stream.write('\n')
destination=R/'attestations/qualification'/revision[:8]/'native-cpu-field-simd-benefit-timeout'
assert not destination.exists();destination.mkdir(parents=True)
sources={str(path.relative_to(R)):path for root in (W,C,D) for path in root.rglob('*') if path.is_file()}
for folder in ('cpp','ispc','generated'):
    for path in (R/'bin/linux-gcc/native-cpu-field-simd'/folder).iterdir():
        assert path.is_file()
        sources[str(path.relative_to(R))]=path
for path in (Path(__file__).resolve(),basis):
    sources[str(path.relative_to(R))]=path
products=[]
for kind in ('cpp','ispc','adapter'):
    products += [R/row['object']['path'] for row in doc(W/(kind+'-compile-result.json'))['objects']]
provider=doc(W/'provider/provider-build-result.json')
products += [R/provider[key] for key in ('library','ICD')]
diagnostic=doc(D/'build-completion.json')
products += [R/diagnostic[key]['path'] for key in ('shim_object','executable')]
assert len(products)==len(set(products))==11
bundles=[archive('source-build-and-trial.zip',sources),archive('actual-products.zip',{str(path.relative_to(R)):path for path in products})]
manifest={'schema':'sirius-native-cpu-field-simd-terminal-v1','qualifying_source_revision':revision,
    'outcome':outcome,'owners':owners,'archives':bundles,'actual_new_product_members':11,
    'source_build_only':False,'full_scientific_or_release_qualification':False,
    'scope':'All9380canonical instructions retained with source dependency replay; actual all3CPP/all3ISPC/adapter/provider builds and exact533function joins PASS, completeTransport38fixture/all187960words and C8/all32/51markers/54joins PASS. Original benefit B1 timed out; no ABBA22/Decimal/all8/frame acceptance or scientific/main/adoption claim.',
    'retained_failed_accounting_owner':'Diagnostic build owner failed only in the old code reader after its single successful shim compile/link; both commands preserved, read-only proof completion passed without rebuilding. No failed owner relabelled.',
    'baseline_and_driver_recovery':'attestations/qualification/603a155e/native-cpu-accumulate-inline-frame-timeout/MANIFEST.json',
    'current_zlib_runtime':doc(C/'runtime-transition.json'),
    'cleanup_scope':'No frame/resource restoration occurred. Failed-candidate products may be retired only after independent archive verification and committed/pushed recovery; active caches and unrelated work remain.'}
with (destination/'MANIFEST.json').open('x') as stream:
    json.dump(manifest,stream,indent=2);stream.write('\n')
print(json.dumps({'bundle':str(destination.relative_to(R)),'owners':len(owners),
    'archives':[{'file':row['file'],'bytes':row['identity']['bytes'],'sha256':row['identity']['sha256'],'members':len(row['members'])} for row in bundles]}))
