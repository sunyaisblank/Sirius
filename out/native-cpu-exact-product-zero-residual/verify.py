"""Reuse the bounded optional compile and original tagged fixture unchanged."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-exact-product-residual/verify.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '4cd9935952f039ab9194b427f01d7fa534d6f3ca6b89b1ca9cc7703cf781ebbf'
source = basis.read_text().replace('native-cpu-exact-product-residual', 'native-cpu-exact-product-zero-residual')
exec(compile(source, str(basis), 'exec'))
