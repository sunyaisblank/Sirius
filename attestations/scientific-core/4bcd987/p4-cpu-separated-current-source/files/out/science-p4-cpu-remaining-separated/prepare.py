"""Seal three separate original CPU controls; never execute them here."""
from pathlib import Path
import ast
import copy
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[1]
PRIOR = ROOT / 'out/science-p4-cpu-remainder'
RUNNER = 'f5d649bf5067d1d68dc9215462f37830122f0540'
PRODUCT = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
LIMITS = {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
CASES = {
    'shakura-sunyaev-live': 'GeodesicTracerTest.LiveDiskTemperatureUsesZeroTorqueShakuraSunyaevProfile',
    'page-thorne-live': 'GeodesicTracerTest.LiveDiskTemperatureUsesFullPageThorneProfile',
    'near-extremal-live': 'GeodesicTracerRedshift.NearExtremalInnerDiskEmissionRemainsFinite',
}
EXCLUDED = 'CpuTraceBoundary.TruncatedPageThorneLiveProfileUsesDeclaredZeroTorqueEdge'
OLD_CASE = 'original-four-non-frame-cpu-controls'


def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def key(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def write(name, value):
    path = WORK / name
    assert not path.exists()
    path.write_text(json.dumps(value, indent=2) + '\n')


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, timeout=30)


assert Path.cwd() == ROOT and {p.name for p in WORK.iterdir()} == {'prepare.py'}
assert git('rev-parse', 'HEAD').decode().strip() == RUNNER and not git('status', '--porcelain').strip()
base = json.loads((PRIOR / 'preparation.json').read_text())
assert base['runner_source_revision'] == RUNNER and base['product_source_revision'] == PRODUCT
assert set(base['expected_names']) == set(CASES.values()) | {EXCLUDED}
for name, pin in base['pins'].items():
    assert identity(ROOT / name) == pin, name
owner = json.loads((PRIOR / OLD_CASE / 'owner.json').read_text())
before = json.loads((PRIOR / OLD_CASE / 'inputs-before.json').read_text())
assert before == json.loads((PRIOR / OLD_CASE / 'inputs-after.json').read_text())
assert owner['status'] == 'terminal' and owner['stop_reason'] == 'diagnostic_timeout' and owner['child_exit'] == -15
assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and not owner['cleanup_errors']
assert not owner['remaining_owned_processes'] and all(row['absent'] for row in owner['observed_births'])
stdout = (PRIOR / OLD_CASE / 'stdout.log').read_text()
for text in ('[       OK ] ' + EXCLUDED + ' (296730 ms)', '[  PASSED  ] 1 test.',
             EXCLUDED + ' ...   Passed  296.75 sec'):
    assert text in stdout
assert not (PRIOR / OLD_CASE / 'ctest.xml').exists() and not (PRIOR / OLD_CASE / 'ctest.log').exists()
original_ids = {str(p.relative_to(PRIOR)): identity(p) for p in PRIOR.rglob('*') if p.is_file()}
pins = {name: pin for name, pin in base['pins'].items() if not (ROOT / name).is_relative_to(PRIOR)}
for name in ('run_owned.py', 'owned_guard.py'):
    expected = {k: base['owner_reuse'][name][k] for k in ('bytes', 'sha256')}
    assert identity(PRIOR / name) == expected
    shutil.copy2(PRIOR / name, WORK / name)
for path in sorted((PRIOR / 'reference').rglob('*')):
    if path.is_file():
        target = WORK / path.relative_to(PRIOR)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        pins[key(target)] = identity(target)
for relative in ('preparation.json', 'actions.json', 'result.json',
                 OLD_CASE + '/owner.json', OLD_CASE + '/stdout.log', OLD_CASE + '/stderr.log',
                 OLD_CASE + '/inputs-before.json', OLD_CASE + '/inputs-after.json',
                 OLD_CASE + '/ctest.log.tmp7d2d4'):
    source, target = PRIOR / relative, WORK / 'prior-interrupted-observation' / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    assert identity(target) == identity(source)
    pins[key(target)] = identity(target)

# Retain the mature runner; only selection, destinations and serial scheduling
# change. Numerical assertions stay in the original executable.
wrapper = (PRIOR / 'run_remainder.py').read_text()
assert wrapper.count('CASE = "' + OLD_CASE + '"') == 1
wrapper = wrapper.replace('One frozen four-case original P4 non-frame remainder action',
                          'Three separate one-shot original P4 CPU actions')
wrapper = wrapper.replace('import hashlib\n', 'import fcntl\nimport hashlib\n')
wrapper = wrapper.replace('CASE = "' + OLD_CASE + '"',
    'CASE_NAMES = ' + repr(CASES) + '\nCASE = sys.argv[1] if len(sys.argv) == 2 else None')
node = next(n for n in ast.parse(wrapper).body if isinstance(n, ast.Assign) and
            any(isinstance(t, ast.Name) and t.id == 'NAMES' for t in n.targets))
lines = wrapper.splitlines(keepends=True)
lines[node.lineno - 1:node.end_lineno] = ['NAMES = {CASE_NAMES.get(CASE)}\n']
wrapper = ''.join(lines)
wrapper = wrapper.replace('(WORK / "result.json")', '(WORK / (CASE + "-result.json"))')
wrapper = wrapper.replace('temporary = WORK / "temporary"', '''require(CASE in CASE_NAMES, "choose one declared original case")
    serial = (WORK / "serial.lock").open("a")
    fcntl.flock(serial, fcntl.LOCK_EX | fcntl.LOCK_NB)
    for other in CASE_NAMES:
        folder = WORK / other
        if other == CASE or not folder.exists():
            continue
        prior_owner = json.loads((folder / "owner.json").read_text())
        require(prior_owner["status"] == "terminal" and not prior_owner["remaining_owned_processes"]
                and not prior_owner["cleanup_errors"] and prior_owner["source_unchanged"]
                and prior_owner["whole_inputs_unchanged"] and prior_owner["observed_births"]
                and all(row["absent"] for row in prior_owner["observed_births"]),
                "another action has not retired cleanly")
    temporary = WORK / ("temporary-" + CASE)''')
wrapper = wrapper.replace('"scope": prepared["scope"], "finite_acceptance_scope": prepared["finite_acceptance_scope"],',
    '"scope": prepared["case_scopes"][CASE], "finite_acceptance_scope": {name: prepared["finite_acceptance_scope"][name] for name in NAMES},')
wrapper = wrapper.replace('set(prepared["expected_names"]) == NAMES and len(NAMES) == 4',
    'set(prepared["expected_names"]) == set(CASE_NAMES.values()) and len(NAMES) == 1')
wrapper = wrapper.replace('set(actions) == {CASE} and actions[CASE]["limits"] == LIMITS',
    'set(actions) == set(CASE_NAMES) and all(action["limits"] == LIMITS for action in actions.values())')
wrapper = wrapper.replace('len(cases) == 4', 'len(cases) == 1')
ast.parse(wrapper, filename='run_remaining.py')
(WORK / 'run_remaining.py').write_text(wrapper)
prepared = copy.deepcopy(base)
prepared['expected_names'] = sorted(CASES.values())
prepared['finite_acceptance_scope'] = {name: base['finite_acceptance_scope'][name] for name in CASES.values()}
prepared.update(prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), executed=False,
    scope='Three separate unchanged original finite CPU P4 controls; serial one-case owners with300s/4GiB/2GiB each. No completed truncated-edge repeat, full frame, device, full P4/science/Mandatory or release claim.',
    ledger_scope={'criterion': 'P4', 'selected_original_names': sorted(CASES.values()), 'selected_count': 3,
                  'completed_truncated_edge_repeated': False},
    case_scopes={case: 'One unchanged original CPU P4 witness: ' + name for case, name in CASES.items()},
    private_tmpdir={'paths': {case: str(WORK / ('temporary-' + case)) for case in CASES},
        'created_only_at_execution': True, 'mode': '0700', 'owner_uid': os.getuid(),
        'preserve_natural_residuals': True},
    serial_execution={'mechanism': 'One nonblocking exclusive flock held through wrapper/owner completion; fail-closed retirement check for any previous action.',
        'lock_path': str(WORK / 'serial.lock'), 'case_order_recommendation': list(CASES),
        'prior_pass_required_for_next': False, 'prior_clean_retirement_required': True})
join = prepared['source_and_registration_join']
join['test_sources'] = ['tests/render/geodesic_tracer_test.cpp']
join['producers'] = {'render': base['source_and_registration_join']['producers']['render']}
join['original_compile_commands'] = [row for row in join['original_compile_commands'] if row['file'].endswith('/tests/render/geodesic_tracer_test.cpp')]
join['original_link_ninja_lines'] = [line for line in join['original_link_ninja_lines'] if line.startswith('build src/sirius/app/sirius_render_tests ')]
join['registrations'] = [row for row in join['registrations'] if row['name'] in CASES.values()]
join['selected_backend_cases'], join['selected_render_cases'] = 0, 3
assert len(join['original_compile_commands']) == len(join['original_link_ninja_lines']) == 1
assert len(join['registrations']) == 3 and all('TIMEOUT' not in row['properties'] for row in join['registrations'])
model = json.loads((ROOT / 'tests/operating_model.json').read_text())
assert set(CASES.values()) <= {x['name'] for x in next(c for c in model['acceptance_criteria'] if c['id'] == 'P4')['evidence']}
for source, saved in base['reference']['original_source_snapshots'].items():
    assert git('show', PRODUCT + ':' + source) == (ROOT / source).read_bytes()
    assert identity(ROOT / source) == {k: saved[k] for k in ('bytes', 'sha256')}
prepared['reference']['cost_scope'] = {
    'full_frame_publications': 0, 'per_case': {
        CASES['shakura-sunyaev-live']: {'maximum_original_Trace_calls': 576, 'max_steps_per_trace': 5000},
        CASES['page-thorne-live']: {'maximum_original_Trace_calls': 576, 'max_steps_per_trace': 5000},
        CASES['near-extremal-live']: {'maximum_original_Trace_calls': 1012, 'max_steps_per_trace': 20000}},
    'scope': 'Original raster searches stop by their original crossing counters; source work bounds are not timing forecasts. Each imported outer300s budget can produce an incomplete observation.'}
prepared['reference']['excluded_completed_component'] = {
    'name': EXCLUDED, 'classification': 'COMPONENT_PASS_OBSERVED_IN_INCOMPLETE_BATCH',
    'source_preparation': key(WORK / 'prior-interrupted-observation/preparation.json'),
    'stdout': key(WORK / 'prior-interrupted-observation' / OLD_CASE / 'stdout.log'),
    'raw_owner': key(WORK / 'prior-interrupted-observation' / OLD_CASE / 'owner.json'),
    'temporary_ctest_log': key(WORK / 'prior-interrupted-observation' / OLD_CASE / 'ctest.log.tmp7d2d4'),
    'producer': base['source_and_registration_join']['producers']['backend'],
    'scope': 'Original first case fully marked GTestOK/PASSED and CTestPassed296.75s; batch remains incomplete. No synthetic XML, standalone passing owner, complete batch or qualification claim.'}
prepared['reference']['gap_boundary'] = 'Only the started-but-interrupted Shakura case and two unstarted original cases are scheduled; completed first component and other passing/reused controls excluded.'
for path in (WORK / 'prepare.py', WORK / 'run_remaining.py', WORK / 'run_owned.py', WORK / 'owned_guard.py'):
    pins[key(path)] = identity(path)
ctest = base['pins']['/usr/bin/ctest']
assert identity(Path('/usr/bin/ctest')) == ctest
actions = {}
for case, name in CASES.items():
    actions[case] = {'revision': RUNNER, 'kind': 'check', 'expected_exit': 0, 'limits': LIMITS,
        'scope': prepared['case_scopes'][case],
        'argv': ['/usr/bin/ctest', '--test-dir', 'bin/linux-gcc', '-V', '--output-on-failure', '--no-tests=error',
            '--output-junit', '{folder}/ctest.xml', '--output-log', '{folder}/ctest.log', '-R', '^' + re.escape(name) + '$'],
        'frozen_inputs': list(pins) + [key(WORK / 'actions.json'), key(WORK / 'preparation.json')]}
write('actions.json', actions)
pins[key(WORK / 'actions.json')] = identity(WORK / 'actions.json')
prepared['pins'] = pins
prepared['owner_reuse']['copied_from_interrupted_root'] = str(PRIOR)
write('preparation.json', prepared)
for name, pin in original_ids.items():
    assert identity(PRIOR / name) == pin
for path in WORK.glob('*.py'):
    ast.parse(path.read_text(), filename=str(path))
print(json.dumps({'prepared_without_execution': True, 'cases': CASES, 'pins': len(pins),
    'owner_guard_byte_equal': all((WORK / n).read_bytes() == (PRIOR / n).read_bytes() for n in ('run_owned.py', 'owned_guard.py')),
    'preparation': identity(WORK / 'preparation.json'), 'runner': identity(WORK / 'run_remaining.py'),
    'commands': [['python3.12', '-B', key(WORK / 'run_remaining.py'), case] for case in CASES]}, indent=2))
