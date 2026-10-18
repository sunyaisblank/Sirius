"""Read completed offline products only; ROOT must admit the owner/input proof first."""
from pathlib import Path
import hashlib
import json
import re
import struct

W = Path(__file__).resolve().parent
R = W.parents[1]
B = R / 'bin/linux-gcc/src/sirius/backend/retained'
C = W / 'candidate'
TARGETS = {'retained_transport.spv': 'kTransportShader',
           'retained_endpoint.spv': 'kEndpointShader'}

def document(path):
    return json.loads(path.read_text())

def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

def pin(path):
    return {'path': str(path.relative_to(R)), **identity(path)}

def require_pin(path, expected):
    assert identity(path) == {key: expected[key] for key in ('bytes', 'sha256')}, str(path)
    if 'resolved_path' in expected:
        assert str(path.resolve()) == expected['resolved_path'], str(path)

def arrays(path):
    result = {}
    pattern = r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};'
    for count, name, body in re.findall(pattern, path.read_text(), re.S):
        words = [int(word) for word in re.findall(r'(\d+)u', body)]
        assert len(words) == int(count) and name not in result
        result[name] = struct.pack('<' + str(len(words)) + 'I', *words)
    assert len(result) == 40
    return result

def main():
    output = W / 'root-compiler-products.json'
    assert not output.exists(), 'Do not overwrite a prior product receipt'
    preflight = document(W / 'root-preflight.json')
    action = document(W / 'actions.json')['candidate']
    owner = document(C / 'owner.json')
    before = document(C / 'inputs-before.json')
    after = document(C / 'inputs-after.json')
    assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['stop_reason'] and not owner['cleanup_errors']
    assert not owner['remaining_owned_processes'] and not owner['provider_mappings']
    assert all(birth['absent'] for birth in owner['observed_births'])
    assert owner['source'] == {'revision': preflight['source_revision'], 'status': ''}
    assert owner['limits'] == action['limits'] == preflight['limits']
    assert owner['argv'] == [arg.replace('{folder}', str(C)) for arg in action['argv']]
    assert before == after and len(before) == preflight['expected_owner_input_rows']
    normalized = {str((R / path).resolve()): row for path, row in before.items()}
    expected_paths = set(preflight['files']) | {str(W / 'root-preflight.json')}
    assert {str(R / path) for path in before} == expected_paths
    # ROOT checks the entire current frozen set; this reader additionally pins all
    # baseline bytes used by its comparisons to the actual precompile input map.
    def frozen(path):
        require_pin(path, normalized[str(path.resolve())])
    baseline_paths = sorted(Path(path) for path in preflight['files']
                            if Path(path).parent == B and path.endswith('.spv'))
    assert len(baseline_paths) == 35
    assert {path.name for path in C.glob('*.spv')} == {path.name for path in baseline_paths}
    products = []
    changed = []
    for baseline in baseline_paths:
        frozen(baseline)
        candidate = C / baseline.name
        old, new = baseline.read_bytes(), candidate.read_bytes()
        assert len(old) % 4 == len(new) % 4 == 0
        assert old[:4] == new[:4] == b'\x03\x02\x23\x07'
        equal = old == new
        if not equal:
            changed.append(baseline.name)
        products.append({'name': baseline.name, 'baseline': identity(baseline),
                         'candidate': identity(candidate), 'byte_exact': equal})
    assert set(changed) == set(TARGETS) and len(changed) == 2
    baseline_header, candidate_header = B / 'retained_kernels.h', C / 'retained_kernels.h'
    frozen(baseline_header)
    old_arrays, new_arrays = arrays(baseline_header), arrays(candidate_header)
    assert set(old_arrays) == set(new_arrays)
    array_changes = sorted(name for name in old_arrays if old_arrays[name] != new_arrays[name])
    assert array_changes == sorted(TARGETS.values())
    for filename, array_name in TARGETS.items():
        assert old_arrays[array_name] == (B / filename).read_bytes()
        assert new_arrays[array_name] == (C / filename).read_bytes()
    programs = sorted(name for name in old_arrays if name.endswith('Program'))
    assert len(programs) == 7
    assert all(old_arrays[name] == new_arrays[name] for name in programs)
    scalar_pattern = r'inline constexpr (std::size_t|bool) (\w+) = ([^;]+);'
    constants = re.findall(scalar_pattern, baseline_header.read_text())
    assert constants and constants == re.findall(scalar_pattern, candidate_header.read_text())
    owner_products = []
    for name in ('retained_owner_kernels.h', 'retained_owner_layout.h', 'retained_owner_layout.slang'):
        frozen(B / name)
        assert (B / name).read_bytes() == (C / name).read_bytes(), name
        owner_products.append({'name': name, 'identity': identity(C / name)})
    result = {'schema': 'sirius-native-fp32-32-lanes-compiler-product-check-v1',
              'passed': True, 'source_revision': preflight['source_revision'],
              'reader': pin(Path(__file__)), 'owner': pin(C / 'owner.json'),
              'before': pin(C / 'inputs-before.json'), 'after': pin(C / 'inputs-after.json'),
              'input_rows': len(before), 'root_whole_current_input_check_required': True,
              'products': products, 'changed_spvs': sorted(changed), 'unchanged_spvs': 33,
              'baseline_header': pin(baseline_header), 'candidate_header': pin(candidate_header),
              'main_arrays': 40, 'changed_arrays': array_changes, 'unchanged_arrays': 38,
              'seven_tables': programs, 'public_constants_exact': constants,
              'owner_products_byte_exact': owner_products,
              'emitted_semantic_comparison_passed': False, 'runtime_released': False,
              'scope': 'Completed compiler/validator products and exact binary exclusions only. '
                       'No emitted ABI/CFG/semantic proof, provider/numerical/control/ABBA/cost/'
                       'adoption/frame/science qualification. ROOT invokes only after owner/input admission.'}
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'products_passed': True, 'changed': sorted(changed),
                      'unchanged': 33, 'runtime_released': False}))

if __name__ == '__main__':
    main()
