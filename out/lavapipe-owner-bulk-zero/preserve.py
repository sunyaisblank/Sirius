"""Reuse terminal preservation with explicit modes and local-build limitations."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/lavapipe-owner-machine-codegen/preserve.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == 'c3c4e0470a673f45edae0bcc2f7fb8bf90e09ba3446c7664c47b30bc5e0e8012'
_source = _basis.read_text()
assert 'lavapipe-owner-machine-codegen-finite-outcome' in _source
_source = _source.replace('lavapipe-owner-machine-codegen-finite-outcome','lavapipe-owner-bulk-zero-finite-outcome')
assert "'machine-codegen.patch'" in _source
_source = _source.replace("'machine-codegen.patch'","'bulk-zero.patch'")
assert "'benefit_attempted':bool(benefit)" in _source
_source = _source.replace("'benefit_attempted':bool(benefit)","'benefit_attempted':any((W/'runtime'/label/'owner.json').exists() for label in ['a1','b1','b2','a2'])")
assert "'sole_ABBA':{'attempted':bool(benefit)" in _source
_source = _source.replace("'sole_ABBA':{'attempted':bool(benefit)","'sole_ABBA':{'attempted':any((W/'runtime'/label/'owner.json').exists() for label in ['a1','b1','b2','a2'])")
assert "'diagnostics':{str(p.relative_to(R)):pin(p) for p in diagnostics}" in _source
_source = _source.replace("'diagnostics':{str(p.relative_to(R)):pin(p) for p in diagnostics}","'diagnostics':{str(p.relative_to(R)):{**pin(p),'mode':stat.S_IMODE(p.stat().st_mode),**({'symlink':os.readlink(p)} if p.is_symlink() else {})} for p in diagnostics}")
assert "'actual_product_paths':paths,'content_addressed_blobs':blobs," in _source
_source = _source.replace("'actual_product_paths':paths,'content_addressed_blobs':blobs,","'actual_product_paths':paths,'content_addressed_blobs':blobs,'sealed_hash_only_reused_local_Mesa_build_inputs':sealed_build,")
assert 'The12s partial stack did not establish the300s timeout cause. No machine-policy selection implies numerical or performance success.' in _source
_source = _source.replace('The12s partial stack did not establish the300s timeout cause. No machine-policy selection implies numerical or performance success.','Exact bulk scratch initialization is only a source contract until emitted retention, resolved publication, full numerical and original benefit gates pass. Original Default machine optimization is retained.')
assert 'External system/toolchain inputs are identified by original before/after manifests and are not all bundled.' in _source
_source = _source.replace('External system/toolchain inputs are identified by original before/after manifests and are not all bundled.','External system/toolchain and separately listed reused local Mesa build inputs are identified by original before/after manifests and are not all bundled. Diagnostic modes/link targets are explicit; ZIP stores target contents, restore relative symlinks after content and mode verification.')
exec(compile(_source,str(_basis),'exec'))
