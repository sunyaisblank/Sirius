"""Complete preservation metadata without altering frozen trial inputs."""
from pathlib import Path
import hashlib
_basis = Path(__file__).resolve().parent / "preserve.py"
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == 'c3c4e0470a673f45edae0bcc2f7fb8bf90e09ba3446c7664c47b30bc5e0e8012'
_source = _basis.read_text()
assert "'benefit_attempted':bool(benefit)" in _source
_source = _source.replace("'benefit_attempted':bool(benefit)", "'benefit_attempted':any((W/'runtime'/label/'owner.json').exists() for label in ['a1','b1','b2','a2'])")
assert "'sole_ABBA':{'attempted':bool(benefit)" in _source
_source = _source.replace("'sole_ABBA':{'attempted':bool(benefit)", "'sole_ABBA':{'attempted':any((W/'runtime'/label/'owner.json').exists() for label in ['a1','b1','b2','a2'])")
assert "'actual_complete_code_review':code," in _source
_source = _source.replace("'actual_complete_code_review':code,", "'actual_complete_code_review':code,'observer_incidents':{'renamed_first_failed_capture':{'historical_path':'ordinary/transport','preserved_path':'ordinary/observer-plt-transport','raw_case_and_argv_unchanged':True,'reason':'Observer resolved provider PLT instead of actual LLVM body; no dispatch/JIT completed, original failed ownership evidence retained.'},'completed_Transport_capture':'Original LLVM C-API creation hook was unused and an unrelated deleted memfd raised a callback warning. Existing completed lookup/emission/greedy/no-fast/no-bridge and exactIR evidence reused; initialDefault is source-composed only. Future E/DP hooks corrected to actual Mesa helper and filtered object mappings; no completed Transport rerun.'},")
assert "'diagnostics':{str(p.relative_to(R)):pin(p) for p in diagnostics}" in _source
_source = _source.replace("'diagnostics':{str(p.relative_to(R)):pin(p) for p in diagnostics}", "'diagnostics':{str(p.relative_to(R)):{**pin(p),'mode':stat.S_IMODE(p.stat().st_mode),**({'symlink':os.readlink(p)} if p.is_symlink() else {})} for p in diagnostics}")
assert "'actual_product_paths':paths,'content_addressed_blobs':blobs," in _source
_source = _source.replace("'actual_product_paths':paths,'content_addressed_blobs':blobs,", "'actual_product_paths':paths,'content_addressed_blobs':blobs,'sealed_hash_only_reused_local_Mesa_build_inputs':sealed_build,")
assert 'External system/toolchain inputs are identified by original before/after manifests and are not all bundled.' in _source
_source = _source.replace('External system/toolchain inputs are identified by original before/after manifests and are not all bundled.', 'External system/toolchain inputs and the separately listed reused local Mesa build inputs are identified by original before/after manifests and are not all bundled. Diagnostic file-link targets/modes are explicit; ZIP stores their authentic target contents. Restore referenced tracked source at its individual owner commit, then diagnostic modes and relative symlinks.')
exec(compile(_source, str(_basis), 'exec'))
