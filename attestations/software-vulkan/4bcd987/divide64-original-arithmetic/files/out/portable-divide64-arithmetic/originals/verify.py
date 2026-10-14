"""Post-owner verification of the unchanged original ordered arithmetic corpus."""
from pathlib import Path
import hashlib
import json
import shlex
import xml.etree.ElementTree as ET
R = Path.cwd().resolve()
W = Path(__file__).resolve().parent
T = R / 'bin/linux-gcc/computed-carry-exit-original-prerequisites/arithmetic'
actions = json.loads((W / 'actions.json').read_text())
def identity(path):
    with Path(path).open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': Path(path).stat().st_size, 'sha256': digest}
for name, pin in json.loads((W.parent / 'source-freeze.json').read_text()).items():
    assert identity(name) == pin, name
for label in actions:
    if label == 'validate-arithmetic':
        continue  # This verifier never reads its own enclosing owner's verdict.
    folder = W / label
    owner = json.loads((folder / 'owner.json').read_text())
    assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and owner['cleanup_errors'] == [] and owner['remaining_owned_processes'] == []
    assert owner['observed_births'] and all(b['absent'] for b in owner['observed_births'])
    assert json.loads((folder / 'inputs-before.json').read_text()) == json.loads((folder / 'inputs-after.json').read_text())
    if label.startswith('execute-gpu-'):
        assert owner['provider_mappings']
        assert (folder / 'actual.bin').read_bytes() == (W / 'expected.bin').read_bytes()
        summary = [json.loads(line) for line in (folder / 'stdout.log').read_text().splitlines() if line.startswith('{')]
        assert len(summary) == 1 and summary[0]['records'] == 174268 and summary[0]['complete_words'] == 697072 and summary[0]['mismatches'] == 0
xml = ET.parse(W / 'execute-host/gtest.xml').getroot()
assert int(xml.get('tests', '-1')) == 1 and not any(xml.findall('.//' + tag) for tag in ['failure', 'error', 'skipped'])
cases = xml.findall('.//testcase')
assert len(cases) == 1 and cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
assert cases[0].get('classname') + '.' + cases[0].get('name') == 'PortableBinary32.MatchesIndependentExactCorpus'
dependencies = shlex.split((T / 'host.d').read_text().replace('\\\n', '').split(':', 1)[1])
candidate_header = W.parent / 'candidate/src/sirius/kernels/portable_binary32.h'
assert candidate_header in {(R / p).resolve() for p in dependencies}
before = json.loads((W / 'compile-host/inputs-before.json').read_text())
for dependency in dependencies:
    path = (R / dependency).resolve()
    key = str(path.relative_to(R)) if path.is_relative_to(R) else str(path)
    assert key in before and identity(path) == {k: before[key][k] for k in ['bytes', 'sha256']}, key
target = W / 'arithmetic-result.json'
assert not target.exists()
target.write_text(json.dumps({'passed': True, 'original_records': 174268, 'complete_output_words_per_provider': 697072,
    'host_dependency_uses_frozen_candidate_header': True, 'full_stock_candidate_outputs_equal_original_expected': True,
    'scope': 'Original unchanged primitive corpus only; focused RTAccumulateTerm phase is separate. Original integrated controls, benefit and scientific qualification remain required.'}, indent=2) + '\n')
print('Original174268 arithmetic corpus accepted for host and both providers.')
