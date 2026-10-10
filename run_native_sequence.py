"""One preregistered serial sequence; stop at the first owner/contract failure."""
from pathlib import Path
import json, os, subprocess, sys, time, signal
from native_contracts import ROOT, WORK, CASES, document, read_case, sha

controller_signals=[];relay_errors=[];active_process=None
def interrupted(signum,frame):
    controller_signals.append(signum)
    if active_process is not None and active_process.poll() is None:
        try:os.kill(active_process.pid,signal.SIGTERM)
        except ProcessLookupError:pass
        except BaseException as error:relay_errors.append(repr(error))
signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)

def owned_call(arguments):
    global active_process
    assert not controller_signals, ('interrupted before spawn', controller_signals)
    call={'arguments':arguments,'cleanup_errors':[],'ownership_verified':False}
    record.setdefault('owned_calls',[]).append(call)
    save()
    active_process=subprocess.Popen(arguments,cwd=ROOT)
    try:
        call['pid']=active_process.pid
        call['start_ticks']=Path('/proc',str(active_process.pid),'stat').read_text().rsplit(')',1)[1].split()[19]
        save()
        if controller_signals and active_process.poll() is None:
            os.kill(active_process.pid,signal.SIGTERM)
        code=active_process.wait()
    except BaseException as error:
        call['cleanup_errors'].append(repr(error))
        if active_process.poll() is None:
            try:active_process.terminate()
            except BaseException as error:call['cleanup_errors'].append(repr(error))
        try:code=active_process.wait(timeout=420)
        except BaseException as error:
            call['cleanup_errors'].append(repr(error));code=active_process.returncode
    call['returncode']=code
    try:current_birth=Path('/proc',str(active_process.pid),'stat').read_text().rsplit(')',1)[1].split()[19]
    except (FileNotFoundError,ProcessLookupError):current_birth=None
    except BaseException as error:
        call['cleanup_errors'].append(repr(error));current_birth=call.get('start_ticks')
    call['current_start_ticks']=current_birth
    call['ownership_verified']=active_process.poll() is not None and call.get('start_ticks') is not None and current_birth!=call['start_ticks']
    call['controller_signals']=list(controller_signals)
    call['relay_errors']=list(relay_errors)
    save()
    if call['ownership_verified']:active_process=None
    # A failed reap retains its live object and exact recorded birth. The sequence
    # reports unverified ownership and cannot start a dependent command.
    return code if call['ownership_verified'] and not call['cleanup_errors'] and not controller_signals and not relay_errors else 1

mode, = sys.argv[1:]
assert mode in ('controls', 'matched')
bindings = document(WORK / 'execution-bindings.json')
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == bindings['candidate_revision']
assert not subprocess.check_output(['git', 'status', '--porcelain'])
if mode == 'controls':
    commands = [('candidate', case) for case in CASES]
else:
    controls = document(WORK / 'native-controls-analysis.json')
    assert controls['pass'] and len(controls['results']) == 2
    for case in CASES:
        read_case('controls', 'candidate', case, bindings)
    assert controls['additive_controls']==[]
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
    result_code = owned_call([sys.executable, str(WORK / 'run_native.py'), mode, label, case])
    if result_code:
        record['failure'] = {'label': label, 'case': case, 'returncode': result_code}
        save(); sys.exit(result_code)
    try:
        observed = read_case(mode, label, case, bindings)
        audit_code=owned_call([sys.executable,str(WORK/'run_terminal_audit.py'),mode,label,case])
        assert audit_code==0
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
assert not controller_signals,('controller signals',controller_signals)
record['completed'] = True
save()
if mode == 'controls':
    analysis = WORK / 'native-controls-analysis.json'
    assert not analysis.exists()
    analysis.write_text(json.dumps({'candidate_revision': bindings['candidate_revision'], 'pass': True, 'results': [q for q in record['results'] if q['case']!='mixed'], 'additive_controls':[q for q in record['results'] if q['case']=='mixed'], 'total_cases':len(record['results']),
                                    'scope': 'Two exact original finite native controls and owner/product/provider/bounded terminal seals. No speed/full-frame/scientific-estate/release claim.'}, indent=2) + '\n')
print(json.dumps({'mode': mode, 'completed': True, 'runs': len(commands)}), flush=True)
