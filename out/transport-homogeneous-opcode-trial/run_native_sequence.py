"""One preregistered serial sequence; stop at the first owner/contract failure."""
from pathlib import Path
import json, os, subprocess, sys, time
from native_contracts import ROOT, WORK, CASES, document, read_case, sha

mode, = sys.argv[1:]
assert mode in ('controls', 'matched')
bindings = document(WORK / 'execution-bindings.json')
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == bindings['candidate_revision']
assert not subprocess.check_output(['git', 'status', '--porcelain'])
if mode == 'controls':
    commands = [('candidate', case) for case in CASES]
else:
    controls = document(WORK / 'native-controls-analysis.json')
    assert controls['pass'] and len(controls['results']) == 9
    for case in CASES:
        read_case('controls', 'candidate', case, bindings)
    assert len(controls['additive_controls'])==1 and controls['additive_controls'][0]['case']=='mixed'
    commands = [(label, 'timestamps') for label in ('a1', 'b1', 'b2', 'a2')]
output = WORK / ('native-' + mode + '-sequence.json')
assert not output.exists()
for label, case in commands:
    assert not (WORK / 'native' / mode / label / case).exists()
record = {'kind': 'Frozen serial original native controls, stop on first failure', 'source_revision': bindings['candidate_revision'],
          'linux_owner_pid': os.getpid(), 'linux_owner_start_ticks': Path('/proc/self/stat').read_text().split(') ', 1)[1].split()[19],
          'observer_sha256': sha(Path(__file__)), 'execution_bindings_sha256': sha(WORK / 'execution-bindings.json'),
          'commands': [{'label': label, 'case': case} for label, case in commands], 'results': [], 'completed': False}
def save():
    output.write_text(json.dumps(record, indent=2) + '\n')
save()
for label, case in commands:
    started = time.monotonic()
    result = subprocess.run([sys.executable, str(WORK / 'run_native.py'), mode, label, case], cwd=ROOT)
    if result.returncode:
        record['failure'] = {'label': label, 'case': case, 'returncode': result.returncode}
        save(); sys.exit(result.returncode)
    try:
        observed = read_case(mode, label, case, bindings)
        audit=subprocess.run([sys.executable,str(WORK/'run_terminal_audit.py'),mode,label,case],cwd=ROOT)
        assert audit.returncode==0
        accepted=document(WORK/('terminal-'+mode+'-'+label+'-'+case+'-accepted.json'))
        assert accepted['pass'] and accepted['record']['artifacts']==observed['artifacts']
        observed['terminal_audit_acceptance']=accepted
        if mode == 'matched':
            reference = record['results'][0] if record['results'] else controls['results'][-1]
            assert observed['timestamps']['readbacks'] == reference['timestamps']['readbacks']
            assert observed['timestamps']['identity'] == reference['timestamps']['identity']
            assert observed['modules'] == reference['modules']
            observed['native_baseline_readback_identity_module_prerequisite_pass'] = True
        record['results'].append(observed)
    except Exception as error:
        record['failure'] = {'label': label, 'case': case, 'reason': str(error), 'elapsed_seconds': time.monotonic() - started}
        save(); raise
    save()
record['completed'] = True
save()
if mode == 'controls':
    analysis = WORK / 'native-controls-analysis.json'
    assert not analysis.exists()
    analysis.write_text(json.dumps({'candidate_revision': bindings['candidate_revision'], 'pass': True, 'results': [q for q in record['results'] if q['case']!='mixed'], 'additive_controls':[q for q in record['results'] if q['case']=='mixed'], 'total_cases':len(record['results']),
                                    'scope': 'Nine exact original plus one additive finite native control and owner/product/provider/bounded terminal seals. No speed/full-frame/scientific-estate/release claim.'}, indent=2) + '\n')
print(json.dumps({'mode': mode, 'completed': True, 'runs': len(commands)}), flush=True)
