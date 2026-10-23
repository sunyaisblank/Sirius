"""Reuse exact all-row controls and original gates with explicit current/reused proof scope."""
from pathlib import Path
import hashlib
_basis=Path.cwd()/'out/lavapipe-intergroup-row-owner/make_runtime_actions.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest()=='cff4ebed1d586144ba8d6c9d0769b7c3db2b9bd28c9d7def998b04ecb63c34a4'
_source=_basis.read_text()
_old="for name in inputs: assert (R / name).is_file(), name"
_new="""inputs |= {'attestations/qualification/337d0bae/lavapipe-owner-bulk-zero-finite-outcome/MANIFEST.json',
           'attestations/qualification/337d0bae/lavapipe-owner-bulk-zero-finite-outcome/diagnostics-and-trial.zip',
           'out/lavapipe-owner-bulk-zero/review-admission.json',
           'out/lavapipe-owner-bulk-zero/source-contract.json',
           'out/lavapipe-owner-bulk-zero/bulk-zero.patch',
           'bin/linux-gcc/lavapipe-owner-bulk-zero/provider/libvulkan_lvp.so'}
for name in inputs: assert (R / name).is_file(), name"""
assert _source.count(_old)==1;_source=_source.replace(_old,_new)
_old='reviewed intergroup SIMD provider'
_new='reviewed complete bulk-zero intergroup SIMD provider, with actual current row-T clear/Default/alpha-equivalent-IR evidence and explicitly source-composed prior E/DP code proof; no current E/DP recapture claim'
assert _source.count(_old)==1;_source=_source.replace(_old,_new)
exec(compile(_source,str(_basis),'exec'))
