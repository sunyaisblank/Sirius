"""Compose exact current Transport code with unchanged E/DP source and disclosed prior code."""
from pathlib import Path
import hashlib,json,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent
S=R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
PRIOR=R/'out/native-original-input-ownership'
def doc(p):return json.loads(p.read_text())
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
def terminal(p):
 o=doc(p)
 assert o['passed'] and o['status']=='terminal' and o['source_unchanged'] and o['whole_inputs_unchanged']
 assert not o['remaining_owned_processes'] and not o['cleanup_errors']
 assert o['observed_births'] and all(x['absent'] for x in o['observed_births'])
 return o
def admit():
 review=doc(W/'review-admission.json');assert review['passed']
 assert pin(W/'bulk-zero.patch')==review['prototype_patch'] and pin(W/'source-contract.json')==review['source_contract']
 for k,v in review['source_files'].items():assert pin(S/k)==v,k
 reuse=doc(W/'reused-row-products.json');assert reuse['three_complete_module_spans_present_once_in_final_ELF']
 for k,v in reuse['restored_actual_products'].items():assert pin(R/k)=={j:v[j] for j in ['bytes','sha256']}
 matched=doc(W/'matched-result.json');exe=R/matched['executable']
 assert matched['passed'] and matched['same_actual_executable_for_both_provider_arms']
 assert pin(exe)['sha256']==matched['executable_sha256']=='7f2d1fa88953b954ab116233458a4cb962470f863fda9bbb3a1cad7d4a7d6a26'
 compiler=doc(W/'compiler-result.json')['modules']
 for kind,v in compiler.items():
  module=R/v['module'];assert pin(module)=={j:v[j] for j in ['bytes','sha256']}
  assert exe.read_bytes().count(module.read_bytes())==1,kind
 provider=doc(W/'providers.json')['providers']['candidate']
 for path,field in [('library','library_identity'),('ICD','ICD_identity')]:assert pin(Path(provider[path]))==provider[field]
 terminal(W/'build/owner.json');built=doc(W/'provider-build-result.json')
 assert built['passed'] and built['library_identity']==provider['library_identity']
 code=doc(W/'ordinary/ordinary-code-review.json');t=code['rows']['Transport']
 assert code['actual_provider']['sha256']==provider['library_identity']['sha256']
 assert all(t[k] for k in ['structural_premise_observed','complete_bulk_clear_observed','remaining_optimized_IR_equivalent','original_Default_machine_pipeline_observed'])
 folder=W/'ordinary/transport';owner=terminal(folder/'owner.json');actual=doc(folder/'capture-verification.json')
 assert actual['passed'] and actual['actual_provider_sha256']==provider['library_identity']['sha256']
 assert actual['actual_complete_module']['sha256']==compiler['Transport']['sha256']
 assert actual['callback_restoration_observed'] and actual['actual_entry']['actual_mxcsr']&0xe040==0
 assert owner['provider_mappings'] and all(provider['library'] in line for m in owner['provider_mappings'] for line in m['maps'])
 properties=doc(PRIOR/'provider-history/candidate-default/float-controls-properties.json')
 creation,=doc(PRIOR/'provider-history/candidate-default/device-features.json')['device_creations']
 observed=doc(folder/'float-controls-properties.json');assert observed and all(x['values']==properties[0]['values'] for x in observed)
 current,=doc(folder/'device-features.json')['device_creations']
 assert current['features']['values']==creation['features']['values']
 assert [x['fma16_32_64'] for x in current['feature_chain']]==[x['fma16_32_64'] for x in creation['feature_chain']]
 for label in ['candidate-default','candidate-poison']:
  proof=doc(PRIOR/'reused/prior-proof/lavapipe-ieee32-validation'/label/'numerical-and-boundary-result.json')
  assert proof['numerical_pass'] and proof['all_observed_callback_states_restored']
 # Only the T hash/2295 shape and cache identifiers changed relative to this
 # authenticated complete E/DP code proof. This is source composition, never
 # a claim that those two stages were recaptured on the current provider.
 D=R/'attestations/qualification/337d0bae/lavapipe-owner-bulk-zero-finite-outcome'
 m=doc(D/'MANIFEST.json');a=D/'diagnostics-and-trial.zip';assert pin(a)==m['archives'][a.name]
 member='out/lavapipe-owner-bulk-zero/ordinary/ordinary-code-review.json'
 with zipfile.ZipFile(a) as z:raw=z.read(member)
 v=m['diagnostics'][member];assert len(raw)==v['bytes'] and hashlib.sha256(raw).hexdigest()==v['sha256']
 prior_code=json.loads(raw)
 for kind in ['Endpoint','DopriPhase']:
  row=prior_code['rows'][kind]
  assert all(row[k] for k in ['structural_premise_observed','complete_bulk_clear_observed','remaining_optimized_IR_equivalent','original_Default_machine_pipeline_observed'])
  assert row['complete_module']['sha256']==compiler[kind]['sha256']
 assert doc(W/'source-contract.json')['basis']['manifest']==str((D/'MANIFEST.json').relative_to(R))
 return True
