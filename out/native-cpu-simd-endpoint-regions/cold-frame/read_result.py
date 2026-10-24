"""Reuse the original F frame reader with L gates and genuine ELF reuse proof."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-dispatch-admission/cold-frame/read_result.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '04b2ff93f78aa150dd4921cbe28aa458b21860d07692bb89a9f44d99b4a73dcf'
source = basis.read_text()
changes = [
    ("current = root / 'out/native-cpu-dispatch-admission'", "current = root / 'out/native-cpu-simd-endpoint-regions'"),
    ("build = doc(work / 'build-result.json')", """historical = root / 'out/native-cpu-dispatch-admission/cold-frame'
assert pin(historical / 'build-result.json')['sha256'] == '558c98ca9463708223e20ba0987f8426f803165757a52c266350d736da4284f3'
build = doc(historical / 'build-result.json')
reuse = doc(work / 'reused-frame.json')
restored = doc(work / 'resource-restore.json')
assert reuse['passed'] and reuse['new_links'] == reuse['new_CPP_compilations'] == 0
assert reuse['historical_build'] == pin(historical / 'build-result.json')
assert reuse['historical_link_recipe'] == pin(historical / 'link-recipe.json')
assert reuse['module_joins'] == build['module_joins']
assert restored['status'] == 'complete' and restored['files'] == build['original_resources_copied_exactly']
assert restored['resource_payload_bytes_copied'] == 0
assert reuse['resource_restore'] == pin(work / 'resource-restore.json')
assert reuse['executable_identity'] == build['executable_identity']
assert build['executable_identity']['sha256'] == 'a767656e1ba8faf8eb08b7ace82f2f076f8c7157d9f336ea2529137daecda1a1'
cache, destination = root / restored['cache'], root / restored['destination']
assert {str(path.relative_to(destination)) for path in destination.rglob('*') if path.is_file()} == set(restored['files'])
for name, value in restored['files'].items():
    assert (cache / name).samefile(destination / name)
    assert pin(cache / name) == pin(destination / name) == value"""),
    ("'build': pin(work / 'build-result.json'),", "'historical_build': pin(historical / 'build-result.json'), 'reuse': pin(work / 'reused-frame.json'), 'resource_restore': pin(work / 'resource-restore.json'),"),
]
for before, after in changes:
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
