"""Separate unchanged numerical checks from this failed acquisition, read-only."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib,json,inspect,subprocess
import verify_console_capture_control as original
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent;output=WORK/'native-capture-console'
read=original.read;sha=original.sha
def seal(p):return {'bytes':p.stat().st_size,'sha256':sha(p)}
plan=read(WORK/'capture-plan.json');ref=read(WORK/'capture-reference.json');o=read(output/'owner.json')
b=read(WORK/'native-capture-console-bridge.json');e=read(WORK/'native-capture-console-emergency.json')
assert not b['accepted'] and b['returncode']==1 and not b['outer_stop'] and not b['errors'] and not b['controller_signals']
assert b['linux_child_birth_absent'] and b['emergency_observer']['linux_birth_absent'] and b['emergency_observer']['returncode']==0
assert o['status']=='failed' and o['exception']=='Capture client did not complete successfully'
assert o['returncode']==0 and o['numerical_control_passed'] and not o['outer_stop'] and not o['cleanup_errors']
assert o['capture_client_returncode']==o['capture_client_terminal_returncode']==1 and not o['capture_trace_sealed']
for k in ('owned_process_absent','test_executable_unchanged','loaded_driver_modules_recorded','source_build_gate_unchanged','capture_ready_observed_before_application','intentional_service_stop','owned_service_absent','owned_capture_client_absent','profiler_tools_unchanged','developer_census_completed','console_detached','console_geometry_verified_after_attachment'):assert o[k] is True,k
assert not o['accessible_remaining_developer_processes'] and e['recorded_births_absent'] and not e['errors']
assert o['outer_guard_seconds']==90 and o['rss_guard_bytes']==4294967296
assert o['source_revision']==plan['native_producer_revision'] and o['live_source_revision']==plan['source_revision'] and o['source_build_gate_sha256']==plan['gate_sha256']
assert all(seal(Path(p) if Path(p).is_absolute() else ROOT/p)==v for p,v in b['input_seals'].items())
# Reuse the exact previously reviewed complete numerical section. No capture
# condition or numerical tolerance is rewritten to accept the failed trace.
reader=WORK/'verify_console_capture_control.py';assert sha(reader)=='e2f815ce53bef7962b7201ca3959efc28f3f865b57d1b4cd1289ed7b31d0eaa8'
source=inspect.getsource(original.verify);start="    exe=ROOT/plan['stage']";end="    trace=output/'original-modern-te.rgp'"
assert source.count(start)==source.count(end)==1
section=source[source.index(start):source.index(end)]
scope=vars(original).copy();scope.update(ROOT=ROOT,WORK=WORK,plan=plan,ref=ref,o=o,output=output)
exec(compile('def unchanged_numerical_checks():\n'+section,'reviewed complete numerical section','exec'),scope)
scope['unchanged_numerical_checks']()
panel=next(v for v in o['scoped_processes'] if v['kind']=='owned-RDP-capture')
screens={}
for n in ('initial','ready','final'):
 s=read(output/(n+'-screen.json'));screens[n]=s;count=s['width']*s['rows_read']
 assert o['required_console_width']<=s['width']<=4096 and 2<=s['height']<=32767 and 0<=s['cursor_y']<s['height']-1 and 0<=s['cursor_x']<s['width']
 assert s['rows_read']==s['cursor_y']+1 and s['requested_count']==s['read_count']==len(s['cells'])==count and 0<count<=2097152
 lines=[s['cells'][i:i+s['width']] for i in range(0,count,s['width'])];assert all(line[-1]==' ' for line in lines)
 assert s['text']==''.join(line.rstrip(' ')+'\r\n' for line in lines)
 allowed={o['owner_pid'],panel['pid']}
 if n=='final':allowed.add(o['test_pid'])
 assert o['owner_pid'] in s['attached_pids'] and len(s['attached_pids'])==len(set(s['attached_pids'])) and set(s['attached_pids'])<=allowed
 if n!='final':assert set(s['attached_pids'])==allowed
assert len({(s['width'],s['height']) for s in screens.values()})==1
assert screens['initial']['cursor_y']<=screens['ready']['cursor_y']<=screens['final']['cursor_y']
for m in ('Profiling (RGP) feature enabled (auto-capture enabled)','Waiting for GPU application to connect...'):assert m in screens['ready']['text']
text=screens['final']['text'];canonical=lambda p:p.replace('\\','/').casefold()
connections=original.re.findall(r'(?m)^Connecting to process: (.+) \(PID: ([0-9]+)\)\r?$',text)
assert len(connections)==2 and all(int(pid)==o['test_pid'] and canonical(path)==canonical(o['test_executable']) for path,pid in connections)
for m in ('[Status] Capturing','[Status] Dumping trace data','15%  0.4 / 2.9 MB','Profiling (RGP) trace failed with result: -2','Capture failed (result: -2)','Successfully set clocks to peak','Successfully restored clocks'):assert m in text,m
assert (output/'panel-screen.txt').read_bytes().decode('utf-8-sig')==text+'\r\n'
assert not list(output.glob('*.rgp'))
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==plan['source_revision']==subprocess.check_output(['git','rev-parse','@{u}'],text=True).strip()
assert not subprocess.check_output(['git','status','--porcelain'],text=True)
record={'numerical_control_passed':True,'capture_acquired':False,'original_bridge_accepted':False,'original_owner_status':'failed','complete_readbacks_exact':32,'selected_original_modules_exact':7,'ordered_available_markers':51,'actual_provider_modules_exact':True,'actual_connected_application_pid':o['test_pid'],'actual_connected_application_path_exact':True,'source_revision':plan['source_revision'],'native_producer_revision':plan['native_producer_revision'],'science_elapsed_seconds':o['science_elapsed_seconds'],'science_sampled_peak_rss_bytes':o['peak_sampled_rss_bytes'],'all_frozen_inputs_independently_rechecked_exact':len(b['input_seals']),'recorded_ownership_closed':True,'automatic_peak_clock_request_and_restore_logged':True,'clock_scope_correction':'Earlier preserved controller/plan wording saying no clocks was inaccurate: RDP automatically reported peak clocks and restoration. No explicit clock or driver experiment configuration was supplied. Actual frequencies and restored hardware state were not independently measured. Profiled timings cannot establish ordinary baseline benefit.','observed_failure':'RDP requested a trace, began dumping, reported 15% / 0.4 of 2.9 MB, then collection errors/disconnection/failed abort and generic result -2; RDP exit1. Cause unestablished.','scope':'Original numerical, stream and provider checks passed under profiler activity; trace acquisition failed. No executed instruction cost, optimization benefit, frame or full qualification.','independent_reviewer':'/root/command_owner_review','evidence':{str(p.relative_to(ROOT)):seal(p) for p in [WORK/'native-capture-console-bridge.json',output/'owner.json',output/'gtest.xml',output/'final-screen.json',WORK/'native-capture-console-emergency.json',reader]}}
p=WORK/'native-capture-console-result.json';assert not p.exists();p.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:record[k] for k in ('numerical_control_passed','capture_acquired','complete_readbacks_exact','selected_original_modules_exact','ordered_available_markers','all_frozen_inputs_independently_rechecked_exact','automatic_peak_clock_request_and_restore_logged')}))
