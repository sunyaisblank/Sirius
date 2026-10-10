"""Invoke one source/product/provider/controller-sealed original native control."""
from pathlib import Path
import json, os, subprocess, sys, time, signal
from native_contracts import ROOT, WORK, CASES, document, sha

mode,label,case=sys.argv[1:]
assert (mode, label, case) == ('controls', 'original', 'science')
producer='baseline' if mode=='baseline' or (mode=='matched' and label.startswith('a')) else 'candidate'
binding_path=WORK/('baseline-native-bindings.json' if mode=='baseline' else 'execution-bindings.json')
bindings = document(binding_path)
head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert head == bindings['candidate_revision'] and not subprocess.check_output(['git', 'status', '--porcelain'])
assert document(WORK / 'candidate-native-apparatus-review.json')['pass']

controller_signals=[]
def interrupted(signum,frame):
    controller_signals.append(signum)
signal.signal(signal.SIGTERM,interrupted)
signal.signal(signal.SIGINT,interrupted)

def verify():
    for relative, record in bindings['files'].items():
        path = Path(relative); path = path if path.is_absolute() else ROOT / path
        assert path.stat().st_size == record['bytes'] and sha(path) == record['sha256'], relative

verify()
assert document(WORK / 'candidate-native-preregistration.json')['readback']['body'] == (WORK / 'candidate-native-preregistration.txt').read_text()
stage = ROOT / bindings['producers'][producer]['stage']
output = WORK / 'native' / mode / label / case
bridge = WORK / ('native-' + mode + '-' + label + '-' + case + '-bridge.json')
assert not output.exists() and not bridge.exists()

def win(path):
    return subprocess.check_output(['wslpath', '-w', str(path.resolve())], text=True).strip()

arguments = ['/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe', '-NoProfile', '-NonInteractive',
             '-ExecutionPolicy', 'Bypass', '-File', win(WORK / 'run_native_control.ps1'),
             '-Stage', win(stage / 'tests/backend/Release'), '-Output', win(output), '-Filter', CASES[case],
             '-Seconds', '300', '-RssBytes', '4294967296', '-TemporaryDirectory', win(WORK / 'temp/native' / mode / label / case),
             '-BuildRoot', win(stage), '-ExpectedGateSha256', bindings['producers'][producer]['gate_sha256'], '-LiveRevision', head]
record = {'kind': 'Finite original/native candidate control or frozen matched comparison bridge with source/product/provider/controller pre/post seals',
          'source_revision': bindings['producers'][producer]['source_revision'], 'live_source_revision': head,
          'source_tree_clean_at_start': True, 'linux_owner_pid': os.getpid(),
          'linux_owner_start_ticks': Path('/proc/self/stat').read_text().split(') ', 1)[1].split()[19],
          'arguments': arguments, 'mode': mode, 'label': label, 'case': case, 'producer': producer,
          'execution_bindings_sha256': sha(binding_path),
          'preregistration_readback_sha256': sha(WORK / 'candidate-native-preregistration.json'),
          'observer_sha256': sha(Path(__file__)), 'started_epoch': time.time(), 'qualification_claimed': False, 'bridge_outer_seconds':480, 'bridge_guard_scope':'Includes Windows owner bootstrap/Add-Type/preflight/test/terminal overhead; original scientific case governors unchanged.',
          'input_bindings': bindings['files'], 'completed': False}
bridge.write_text(json.dumps(record, indent=2) + '\n')
cleanup_args = [arguments[0],'-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(WORK/'emergency_native_cleanup.ps1'),
    '-TemporaryDirectory',win(WORK/'temp/native'/mode/label/case),'-Output',win(output),'-Receipt',win(WORK/('native-'+mode+'-'+label+'-'+case+'-emergency-cleanup.json'))]
started = time.monotonic()
bridge_error, cleanup_error, outer_stop = None, None, False
with bridge.with_suffix('.stdout').open('wb') as out, bridge.with_suffix('.stderr').open('wb') as err:
    assert not controller_signals, ('interrupted before spawn',controller_signals)
    child = subprocess.Popen(arguments, stdout=out, stderr=err, start_new_session=True)
    try:
        record['linux_bridge_child_pid'] = child.pid
        record['linux_bridge_child_start_ticks'] = Path('/proc',str(child.pid),'stat').read_text().rsplit(')',1)[1].split()[19]
        bridge.write_text(json.dumps(record,indent=2)+'\n')
        while child.poll() is None:
            if controller_signals:
                outer_stop=True
                bridge_error='controller_signal'
                break
            remaining=480-(time.monotonic()-started)
            if remaining<=0:
                raise subprocess.TimeoutExpired(arguments,480)
            try:child.wait(timeout=min(0.25,remaining))
            except subprocess.TimeoutExpired:pass
        code=child.returncode
    except BaseException as error:
        bridge_error = repr(error)
        outer_stop = True
    finally:
        if controller_signals:
            outer_stop=True
            bridge_error=bridge_error or 'controller_signal'
        if outer_stop or child.poll() is None:
            # The helper never starts/compiles scientific work. It uses exact early
            # Windows births and only accessible task-bound descendant evidence.
            try:
                with (WORK/('native-'+mode+'-'+label+'-'+case+'-emergency-cleanup.stdout')).open('wb') as cleanup_out, (WORK/('native-'+mode+'-'+label+'-'+case+'-emergency-cleanup.stderr')).open('wb') as cleanup_err:
                    cleanup_result = subprocess.run(cleanup_args,stdout=cleanup_out,stderr=cleanup_err,timeout=60)
                record['emergency_cleanup_returncode'] = cleanup_result.returncode
                assert cleanup_result.returncode == 0
                cleanup_document = document(WORK/('native-'+mode+'-'+label+'-'+case+'-emergency-cleanup.json'))
                assert cleanup_document['recorded_births_absent'] and not cleanup_document['errors']
            except BaseException as error:
                cleanup_error = repr(error)
        try:
            code = child.wait(timeout=15)
        except BaseException as error:
            cleanup_error = cleanup_error or repr(error)
            try:
                os.killpg(child.pid,signal.SIGKILL)
            except ProcessLookupError: pass
            except BaseException as signal_error: cleanup_error = repr(signal_error)
            try: code = child.wait(timeout=10)
            except BaseException as last_error:
                cleanup_error = repr(last_error)
                code = child.returncode
        record.update(bridge_error=bridge_error,cleanup_error=cleanup_error,outer_stop=outer_stop)
        saved_birth = record.get('linux_bridge_child_start_ticks')
        try: current_birth = Path('/proc',str(child.pid),'stat').read_text().rsplit(')',1)[1].split()[19]
        except (FileNotFoundError,ProcessLookupError): current_birth = None
        except BaseException as birth_error:
            current_birth = saved_birth
            cleanup_error = repr(birth_error)
            record['cleanup_error'] = cleanup_error
        record['linux_bridge_child_birth_absent'] = saved_birth is not None and current_birth != saved_birth
if controller_signals:
    outer_stop=True
    bridge_error=bridge_error or 'controller_signal'
    record.update(outer_stop=outer_stop,bridge_error=bridge_error)
record['controller_signals']=list(controller_signals)
record.update(returncode=code, elapsed_seconds=time.monotonic() - started, completed=True)
try:
    verify()
    record['input_seals_current_and_exact_at_end'] = True
except Exception as error:
    record['input_seals_current_and_exact_at_end'] = False
    record['post_seal_error'] = str(error)
record['live_source_revision_at_end'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
record['source_tree_clean_at_end'] = not subprocess.check_output(['git', 'status', '--porcelain'])
record['execution_bindings_unchanged_at_end'] = sha(binding_path) == record['execution_bindings_sha256']
record['controller_signals']=list(controller_signals)
if controller_signals:
    record['outer_stop']=True
    record['bridge_error']=record['bridge_error'] or 'controller_signal'
record['accepted_owner_execution'] = (not controller_signals and code == 0 and not outer_stop and bridge_error is None and cleanup_error is None and record['linux_bridge_child_birth_absent'] and record['input_seals_current_and_exact_at_end'] and record['source_tree_clean_at_end'] and record['live_source_revision_at_end'] == head and record['execution_bindings_unchanged_at_end'])
bridge.write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({k: record[k] for k in ('case', 'label', 'source_revision', 'returncode', 'elapsed_seconds', 'accepted_owner_execution')}), flush=True)
sys.exit(0 if record['accepted_owner_execution'] else 1)
