"""Seal current records for the unchanged native timestamp/ABBA readers; no launch."""
from pathlib import Path
import copy
import hashlib
import json
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'out/native-fp32-packed-scratch-runtime-source'
PRODUCER = '45d91ac36d1135d18cba8ebb1923339723fde5df'
RUNNER = '5b2df58bb9f24e7db5f7fbc0c38bdd832a58900e'


def identity(path):
    digest = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
            size += len(block)
    return {'bytes': size, 'sha256': digest.hexdigest()}


def document(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def main():
    phase, = sys.argv[1:]
    assert phase in ('timestamp', 'abba')
    work = ROOT / ('out/native-fp32-packed-soa-runtime-' + phase)
    output = ROOT / ('out/native-fp32-packed-soa-' +
                     ('timestamp-evidence' if phase == 'timestamp' else 'ABBA-evidence'))
    assert not output.exists()
    bindings = document(work / 'execution-bindings.json')
    plan = document(work / 'plan.json')
    assert bindings['source_revision'] == plan['source_revision'] == PRODUCER
    assert bindings['runner_source_revision'] == plan['runner_source_revision'] == RUNNER
    assert bindings['control_profile'] == plan['control_profile'] == phase
    labels = ['control-8'] if phase == 'timestamp' else ['a1', 'b1', 'b2', 'a2']
    assert plan['order'] == labels
    files = {}

    def add(path, expected=None):
        path = Path(path)
        assert path.is_file()
        key = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
        current = identity(path)
        if expected is not None:
            assert current == {k: expected[k] for k in ('bytes', 'sha256')}, key
        assert key not in files or files[key] == current
        files[key] = current

    for key, expected in bindings['files'].items():
        path = Path(key)
        add(path if path.is_absolute() else ROOT / path, expected)
    for directory in (work, ROOT / 'renders' / work.name):
        for path in sorted(directory.rglob('*')):
            if path.is_file():
                assert not path.is_symlink()
                add(path)
    for name in ('native_evidence.py', 'numeric_gates.py', 'compare_native.py',
                 'read_controls_native.py', 'expected-readback-layout.json'):
        add(SOURCE / name)
    add(Path(__file__).resolve())
    authority = ROOT / 'out/native-fp32-packed-soa-pe-audit/matched-array-authority-45d91ac3-reader.json'
    correction = ROOT / 'out/native-fp32-packed-soa-pe-audit/root-reader-authority-key-correction.json'
    add(authority)
    add(correction)
    av = document(authority)
    assert av['passed'] is True and av['producer_source_revision'] == PRODUCER
    assert all(set(a['fma_available']) == {'transport', 'endpoint'} for a in av['arms'].values())
    cases = {}
    loaders = {}
    hardware = None
    for label in labels:
        accepted = work / ('terminal-' + label + '-accepted.json')
        assert document(accepted)['pass'] is True
        arm = 'baseline' if label in ('a1', 'a2') else 'candidate'
        cases[label] = {'work': str(work.relative_to(ROOT)), 'arm': arm,
                        'profile': phase, 'terminal_accepted': identity(accepted)}
        modules = document(work / 'native' / label / 'driver-modules.json')
        loader, = [m for m in modules if m['name'].casefold() == 'vulkan-1.dll']
        assert arm not in loaders or loaders[arm] == loader
        loaders[arm] = loader
        props = {}
        for node in ET.parse(work / 'native' / label / 'gtest.xml').findall('.//property'):
            key = node.get('name')
            assert key not in props
            props[key] = node.get('value')
        observed = {k: props[k] for k in ('selected_device', 'selected_driver',
                    'retained_scalar_route', 'fma_fp32_enabled',
                    'preserves_fp32_signed_zero_inf_nan', 'timestamp_valid_bits', 'timestamp_period_ns')}
        assert hardware is None or hardware == observed
        hardware = observed
    if phase == 'timestamp':
        baseline = copy.deepcopy(loaders['candidate'])
        baseline['path'] = baseline['path'].replace('\\candidate\\', '\\baseline\\')
        stage = ROOT / bindings['stage']
        assert identity(stage / 'payload/baseline/vulkan-1.dll') == identity(stage / 'payload/candidate/vulkan-1.dll')
        loaders['baseline'] = baseline
    # Expected same-provider capability declarations, not evidence of wide execution.
    hardware.update(supports_fp64=True, rounds_fp64_to_nearest=True)
    value = {'schema': 'sirius-native-packed-evidence-bindings-v1', 'admitted': True,
             'phase': phase, 'producer_source_revision': PRODUCER, 'runner_source_revision': RUNNER,
             'observed_instance_api_minimum_improvement_percent': 0, 'files': files,
             'array_authority': {'path': str(authority.relative_to(ROOT)), **identity(authority)},
             'native_identity': hardware,
             'native_identity_scope': 'Current timestamp XML observations. FP64/RTE64 are expected capability declarations for the same bound AMD/SDK provider; original Factored wide control remains required before all8/adoption.',
             'loader_modules': loaders,
             'loader_module_scope': ('Candidate observed in this child. Baseline path is expected for the identical other staged SDK DLL and must be observed in ABBA.' if phase == 'timestamp' else 'Both arms actually observed in the current ABBA children.'),
             'cases': cases}
    output.mkdir()
    context = output / 'reader-bindings.json'
    context.write_text(json.dumps(value, indent=2) + '\n')
    admission = {'schema': 'sirius-native-packed-evidence-root-admission-v1', 'passed': True,
                 'reader_bindings': {'path': str(context.relative_to(ROOT)), **identity(context)},
                 'scope': 'native matched retained controls and sole original22 ABBA only'}
    (output / 'root-admission.json').write_text(json.dumps(admission, indent=2) + '\n')
    print(json.dumps({'phase': phase, 'files': len(files), 'output': str(output.relative_to(ROOT)),
                      'workload_executed': False}))


if __name__ == '__main__':
    assert not sys.flags.optimize
    main()
