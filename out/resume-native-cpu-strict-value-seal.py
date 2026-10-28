"""Finish the interrupted archive; reuse its complete verified source ZIP."""
from pathlib import Path
import hashlib
import json
import textwrap

basis = Path.cwd() / 'out/native-cpu-strict-ir-linking/seal_value.py'
source = basis.read_text()
assert hashlib.sha256(source.encode()).hexdigest() == 'cadb1a68786c447400fe3d3c6990d54bf02aa2d2805205c19717fe8630cd4261'
replacements = {
    'assert not A.exists()': "assert A.is_dir() and {p.name for p in A.iterdir()} == {'source-and-receipts.zip'}",
    "('dp-value.bc', 'dp-value.joins.bc'": "('dp-value.joins.bc',",
    'A.mkdir(parents=True)': '# Existing source ZIP is immutable and reverified below.',
}
# Remove only the nonexistent failed product, never synthesize its bytes.
replacements["('dp-value.bc', 'dp-value.joins.bc'"] = "('dp-value.joins.bc'"
for old, new in replacements.items():
    assert source.count(old) == 1, old
    source = source.replace(old, new, 1)
start = source.index("    with zipfile.ZipFile(target, 'x',")
end = source.index('    with zipfile.ZipFile(target) as zipped:', start)
block = source[start:end]
source = source[:start] + '    if not target.exists():\n' + textwrap.indent(block, '    ') + source[end:]
scope = {'__file__': str(basis), '__name__': 'resumed_value_seal'}
exec(compile(source, str(basis), 'exec'), scope)
manifest_path = scope['A'] / 'MANIFEST.json'
manifest = json.loads(manifest_path.read_text())
manifest['packing_continuation'] = dict(
    reason='First packing attempt completed/verified source ZIP, then encountered absent dp-value.bc from the deliberately refused transform. No failed compiler product existed; its owner/stderr remain in source ZIP.',
    original_producer_sha256='cadb1a68786c447400fe3d3c6990d54bf02aa2d2805205c19717fe8630cd4261',
    verified_source_ZIP_reused=True, compiler_or_workload_reruns=0)
manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
