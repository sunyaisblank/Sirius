"""Validate the one completed original trial and retain its fixed verdict."""
from pathlib import Path
import hashlib
import json
import subprocess
import xml.etree.ElementTree as ET

root = Path.cwd().resolve()
work = Path(__file__).resolve().parent
trial = work / 'benefit'
load = lambda path: json.loads(path.read_text())


def identity(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


comparison = load(trial / 'software-comparison.json')
actions = load(trial / 'actions.json')
freeze = load(trial / 'retention-gate-frozen.json')
assert load(trial / 'software-sequence.json')['completed']
assert freeze['sequence'] == ['a1','b1','b2','a2'] and freeze['within_slot_repeats'] == 3
assert freeze['complete_call_minimum_improvement_percent'] == 5
assert freeze['every_inclusive_stage_minimum_improvement_percent'] == -5
assert freeze['observed_instance_api_minimum_improvement_percent'] == 0
assert freeze['all20_gates_required']
assert comparison['readback_identity_module_and_completion_prerequisites_pass']
assert len(comparison['checks']) == 20 and len(comparison['accounting_checks']) == 2
assert sum(c['scope'] == 'complete_call_ms' for c in comparison['checks']) == 4
assert comparison['pass'] == all(c['pass'] for c in comparison['checks'] + comparison['accounting_checks'])
assert load(work / 'producer-reuse.json')['passed']
assert load(work / 'control-acceptance.json')['passed']
assert load(work / 'endpoint-route-acceptance.json')['passed']
assert not (work / 'trial-acceptance.json').exists() and not (work / 'disposition.json').exists()
rows, verified = [], set()
for label in freeze['sequence']:
    folder = trial / label
    owner = load(folder / 'owner.json')
    assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
    assert not owner['stop_reason'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['provider_mappings'] and owner['observed_births']
    assert all(item['absent'] for item in owner['observed_births'])
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['source'] == {'revision': freeze['candidate_revision'], 'status': ''}
    assert owner['argv'] == actions[label]['argv']
    before = load(folder / 'inputs-before.json')
    assert before == load(folder / 'inputs-after.json')
    for name, pin in before.items():
        assert identity(root / name) == {key: pin[key] for key in ['bytes','sha256']}, name
        verified.add(name)
    xml = ET.parse(folder / 'gtest.xml').getroot()
    assert int(xml.get('tests', '-1')) == 1
    assert not any(xml.findall('.//' + tag) for tag in ['failure','error','skipped'])
    for node in [xml, *xml.findall('.//testsuite')]:
        assert all(int(node.get(key, '0')) == 0 for key in ['failures','errors','skipped','disabled'])
    cases = xml.findall('.//testcase')
    assert len(cases) == 1 and cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
    rows.append({'label': label, 'wall_seconds': owner['wall_seconds'],
                 'peak_sampled_rss_kib': owner['peak_sampled_rss_kib'],
                 'owner': identity(folder / 'owner.json'), 'xml': identity(folder / 'gtest.xml')})
protected = load(work / 'protected-before.json')
for name, pin in protected.items():
    assert identity(root / name) == pin, name
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip() == freeze['candidate_revision']
assert not subprocess.check_output(['git','status','--porcelain'])
counts = {'interval_stage_passed': sum(c['pass'] for c in comparison['checks']),
          'interval_stage_required': 20,
          'inclusive_API_passed': sum(c['pass'] for c in comparison['accounting_checks']),
          'inclusive_API_required': 2}
acceptance = {'passed': True, 'original_trial_complete': True,
              'source_revision': freeze['candidate_revision'], 'slots': rows,
              'declared_current_inputs_verified': len(verified),
              'protected_payloads_unchanged': len(protected), **counts,
              'performance_acceptance': comparison['pass'],
              'scope': 'Original comparator and complete owner/XML/input/source closure. A rejected fixed benefit verdict is a valid experiment outcome, not full qualification.'}
disposition = {'decision': 'PROCEED_TO_PRODUCTION_IMPLEMENTATION' if comparison['pass'] else 'REJECT_OPTIONAL64_DIVIDE_COST',
               'issue': 44, 'source_revision': freeze['candidate_revision'],
               'original_benefit_trial_completed': True, 'original_correctness_prerequisites_passed': True,
               'production_adopted': False, **counts,
               'complete_call_checks': [c for c in comparison['checks'] if c['scope'] == 'complete_call_ms'],
               'inclusive_API_checks': comparison['accounting_checks'],
               'failed_stage_checks': [c for c in comparison['checks'] if not c['pass'] and c['scope'] != 'complete_call_ms'],
               'policy': 'One unchanged original A1/B1/B2/A2; no retries, pooling, cost subtraction or weakened gates. Source/compiler/arithmetic/eight-control prerequisites remain finite evidence. Full owned process wall retained; original device fixture SetUp remains outside original measured windows. Rendering/science/native/platform/release goal remains active.'}
for name, value in [('trial-acceptance.json', acceptance), ('disposition.json', disposition)]:
    (work / name).write_text(json.dumps(value, indent=2) + '\n')
print(json.dumps({'decision': disposition['decision'], **counts,
                  'verified_inputs': len(verified), 'wall_seconds': sum(r['wall_seconds'] for r in rows)}))
