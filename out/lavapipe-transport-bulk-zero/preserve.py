"""Preserve the terminal Transport trial with exact source reconstruction and proof reuse."""
from pathlib import Path
import hashlib
_adapter=Path.cwd()/'out/lavapipe-owner-bulk-zero/preserve_complete.py'
assert hashlib.sha256(_adapter.read_bytes()).hexdigest()=='f27e336a492eb158c0e3ff7bb587e099cae8c4191ea655e98debe65db5da29c8'
_wrapper=_adapter.read_text();_tail="exec(compile(_source, str(_basis), 'exec'))\n"
assert _wrapper.endswith(_tail)
_old="_adapter = Path(__file__).resolve().with_name('preserve.py')"
assert _wrapper.count(_old)==1
_wrapper=_wrapper.replace(_old,"_adapter = Path.cwd()/'out/lavapipe-owner-bulk-zero/preserve.py'")
exec(compile(_wrapper[:-len(_tail)],str(_adapter),'exec'))
for _old,_new in {
 'lavapipe-owner-bulk-zero-finite-outcome':'lavapipe-transport-bulk-zero-finite-outcome',
 'f958a60563b64c79d3f2fb6da2815bde0effdbbe49f121034cfd212ab17b27ff':'7f2d1fa88953b954ab116233458a4cb962470f863fda9bbb3a1cad7d4a7d6a26',
 "'actual_mixed_host_ELF_reused'":"'actual_all_row_host_ELF_reused'",
 "'prior_bundle':'attestations/qualification/263abe74/lavapipe-intergroup-row-owner-finite-outcome/MANIFEST.json'":"'prior_bundle':'attestations/qualification/337d0bae/lavapipe-owner-bulk-zero-finite-outcome/MANIFEST.json'",
 "'prior_manifest_identity':pin(R/'attestations/qualification/263abe74/lavapipe-intergroup-row-owner-finite-outcome/MANIFEST.json')":"'prior_manifest_identity':pin(R/'attestations/qualification/337d0bae/lavapipe-owner-bulk-zero-finite-outcome/MANIFEST.json')",
 "'prior_preserved_at_commit':'206ca0d4'":"'prior_preserved_at_commit':'48adb23d'",
 'Authenticate prior bundle; reconstruct authentic Mesa25.2.8 + exact5 IEEE32 prerequisite + original263 intergroup prototype, then apply bulk-zero.patch with zero fuzz. Source contract and ROOT admission bind the4 changed files and11100-file source identity. The rejected CodeGenNone patch is absent; original Default machine optimization is restored. Configured original build inputs are sealed in the build owner manifests; actual staged provider is a content-addressed product.':'Authenticate prior337 bundle at48adb23d and its reconstruction chain: authentic Mesa25.2.8 + exact5 IEEE32 + original263 intergroup + the prior E/DP bulk-zero.patch. Then apply this bundle’s distinct Transport bulk-zero.patch with zero fuzz. The4-file delta and full11100-file source identity are authenticated. CodeGenNone remains absent and original Default machine optimization remains. Current T code is actually captured; E/DP source/prior-code reuse is disclosed. Configured build inputs remain sealed/hash-only as declared; staged provider is a content-addressed product.',
}.items():
 assert _source.count(_old)==1,(_old,_source.count(_old))
 _source=_source.replace(_old,_new)
_old='for directory,commit in priors:'
assert _source.count(_old)==1
_source=_source.replace(_old,"priors.append(('337d0bae/lavapipe-owner-bulk-zero-finite-outcome','48adb23d'))\n"+_old)
_old='Exact bulk scratch initialization is only a source contract until emitted retention, resolved publication, full numerical and original benefit gates pass. Original Default machine optimization is retained.'
_new='Current exact row-T complete clear, resolved publication, original Default engine/greedy and full remaining-IR alpha equivalence are observed. E/DP code proof is explicitly historical and source-composed; no fresh E/DP capture is claimed. Numerical/benefit/frame admission is determined only by the original gates below.'
assert _source.count(_old)==1;_source=_source.replace(_old,_new)
exec(compile(_source,str(_basis),'exec'))
