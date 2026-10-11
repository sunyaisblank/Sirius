"""Read complete original control results; no device or tool execution."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, math, re, xml.etree.ElementTree as ET
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fields(s):
    out={}
    for part in s.split(';'):
        key,sep,value=part.partition('=');assert sep and key not in out;out[key]=value
    return out
def verify(output):
    plan=read(WORK/'capture-interval-plan.json');ref=read(WORK/'capture-reference.json');o=read(output/'owner.json')
    assert o['status']=='completed' and o['returncode']==0 and not o['outer_stop'] and not o['cleanup_errors']
    for k in ('owned_process_absent','numerical_control_passed','test_executable_unchanged','loaded_driver_modules_recorded','source_build_gate_unchanged','capture_ready_observed_before_application','capture_trace_sealed','intentional_service_stop','owned_service_absent','owned_capture_client_absent','profiler_tools_unchanged','developer_census_completed','console_detached','console_geometry_verified_after_attachment'):assert o[k] is True,k
    assert not o['accessible_remaining_developer_processes'] and o['capture_client_returncode']==0 and o['capture_client_terminal_returncode']==0
    assert type(o['test_pid']) is int and o['test_pid']>0 and o['test_start_time_utc']
    assert o['outer_guard_seconds']==90 and o['rss_guard_bytes']==4294967296
    assert o['source_revision']==plan['native_producer_revision'] and o['live_source_revision']==plan['source_revision'] and o['source_build_gate_sha256']==plan['gate_sha256']
    exe=ROOT/plan['stage']/'tests/backend/Release/sirius_backend_tests.exe'
    assert o['test_executable_sha256']==sha(exe)
    assert len(o['scoped_processes'])==2
    assert {v['kind'] for v in o['scoped_processes']}=={'owned-RDS','owned-RDP-capture'}
    for v in o['scoped_processes']:
        assert type(v['pid']) is int and v['pid']>0 and v['birth']
        n='RadeonDeveloperServiceCLI.exe' if v['kind']=='owned-RDS' else 'RadeonDeveloperPanelCLI.exe'
        assert v['executable'].replace('\\','/').endswith('/'+n) and v['sha256']==sha(ROOT/plan['tools']/n)
    for key,value in o['environment'].items():
        if key not in ('TEMP','TMP'):assert value in ('',None),(key,value)
    xml=output/'gtest.xml';assert 0<xml.stat().st_size<16*1024*1024
    doc=ET.fromstring(xml.read_bytes());assert doc.get('tests')=='1'
    for node in [doc,*doc.findall('.//testsuite')]:
        for k in ('failures','errors','skipped','disabled'):assert int(node.get(k,'0'))==0,k
    cases=doc.findall('.//testcase');assert len(cases)==1
    assert cases[0].get('classname')+'.'+cases[0].get('name')==plan['filter']
    assert cases[0].get('status')=='run' and cases[0].get('result')=='completed'
    assert not any(doc.findall('.//'+k) for k in ('failure','error','skipped'))
    props={}
    for v in doc.findall('.//property'):
        assert v.get('name') not in props;props[v.get('name')]=v.get('value')
    for k,v in ref['provider_properties'].items():assert props[k]==v,(k,props.get(k),v)
    assert {k:v for k,v in props.items() if k.startswith('selected_module_')}==ref['selected_modules']
    assert {k:v for k,v in props.items() if '_readback_' in k}==ref['readbacks']
    assert props['device_observations']=='51'
    assert {k for k in props if k.startswith('device_observation_')}=={'device_observation_'+str(i) for i in range(51)}
    for i,old in enumerate(ref['ordered_marker_labels']):
        m=fields(props['device_observation_'+str(i)])
        assert set(m)==set(old) and m['phase']==old['phase'] and m['stage']==old['stage']
        assert m['available0']==m['available1']=='1'
        start,end=int(m['begin']),int(m['end']);assert 0<=start<2**64 and 0<=end<2**64
        for key in ('device_ms','host_submit_ms','host_wait_ms','host_completion_ms','pipeline_ms','setup_ms','cleanup_ms','total_ms'):
            v=float(m[key]);assert math.isfinite(v) and v>=0;m[key]=v
        assert math.isclose(m['device_ms'],((end-start)&(2**64-1))*10/1e6,rel_tol=1e-9,abs_tol=1e-9)
        assert math.isclose(m['host_completion_ms'],m['host_submit_ms']+m['host_wait_ms'],rel_tol=1e-9,abs_tol=1e-9)
        assert math.isclose(m['total_ms'],m['pipeline_ms']+m['setup_ms']+m['host_completion_ms']+m['cleanup_ms'],rel_tol=1e-9,abs_tol=1e-9)
    mods=read(output/'driver-modules.json')
    normalize=lambda v:{k:(v[k].casefold() if k in ('name','path') else v[k]) for k in ('name','path','bytes','sha256','version')}
    assert sorted(map(normalize,mods),key=lambda x:x['name'])==sorted(map(normalize,ref['driver_modules']),key=lambda x:x['name'])
    trace=output/'original-legacy-interval.rgp';assert 0<trace.stat().st_size<=268435456
    assert trace.stat().st_size==o['trace']['bytes'] and sha(trace)==o['trace']['sha256']
    assert o['trace']['path'].replace('\\','/').endswith('/native-capture-interval-console/original-legacy-interval.rgp')
    log=output/'panel-screen.txt';assert 0<log.stat().st_size<16*1024*1024
    screens={}
    for name in ('initial','ready','final'):
        s=read(output/(name+'-screen.json'));screens[name]=s
        assert o['required_console_width']<=s['width']<=4096 and 2<=s['height']<=32767
        assert 0<=s['cursor_x']<s['width'] and 0<=s['cursor_y']<s['height']-1 and s['rows_read']==s['cursor_y']+1
        count=s['width']*s['rows_read'];assert 0<count<=2097152
        assert s['requested_count']==s['read_count']==count==len(s['cells'])
        lines=[s['cells'][i:i+s['width']] for i in range(0,count,s['width'])]
        assert all(line[-1]==' ' for line in lines)
        assert s['text']==''.join(line.rstrip(' ')+'\r\n' for line in lines)
        assert o['owner_pid'] in s['attached_pids'] and len(set(s['attached_pids']))==len(s['attached_pids'])
        allowed={o['owner_pid'],next(v['pid'] for v in o['scoped_processes'] if v['kind']=='owned-RDP-capture')}
        if name=='final':allowed.add(o['test_pid'])
        assert set(s['attached_pids'])<=allowed
        if name!='final':assert set(s['attached_pids'])==allowed
    assert len({(s['width'],s['height']) for s in screens.values()})==1
    assert screens['initial']['cursor_y']<=screens['ready']['cursor_y']<=screens['final']['cursor_y']
    assert 'Profiling (RGP) feature enabled (auto-capture enabled)' in screens['ready']['text']
    assert 'Waiting for GPU application to connect...' in screens['ready']['text']
    text=log.read_bytes().decode('utf-8-sig')
    assert text==screens['final']['text']+'\r\n'
    text=screens['final']['text']
    connections=re.findall(r'(?m)^Connecting to process: (.+) \(PID: ([0-9]+)\)\r?$',text)
    assert connections,'No actual profiler application identity recorded'
    canonical=lambda p:p.replace('\\','/').casefold()
    assert all(int(pid)==o['test_pid'] and canonical(path)==canonical(o['test_executable']) for path,pid in connections),connections
    completions=re.findall(r'(?m)^Profiling \(RGP\) trace completed successfully! Size: ([0-9]+) bytes\r?$',text)
    assert len(completions)==1 and int(completions[0])==trace.stat().st_size,completions
    saved=re.findall(r'(?m)^Trace saved to: (.+)\r?$',text)
    assert len(saved)==1 and canonical(saved[0].rstrip('\r'))==canonical(o['trace']['path']),saved
    assert 'Waiting for GPU application to connect...' in text
    return {'original_testcase_passed':True,'complete_readbacks_exact':32,'selected_original_modules_exact':7,'ordered_available_markers':51,'actual_provider_modules_exact':True,'actual_connected_application_pid':o['test_pid'],'actual_connected_application_path_exact':True,'trace_completion_and_save_logged':True,'trace_bytes':trace.stat().st_size,'trace_sha256':sha(trace),'scope':'Original numerical/stream/provider control plus whole trace acquisition. Trace dispatch/pipeline identity and sample interpretation remain a separate analysis prerequisite; no performance benefit, frame or full qualification.'}
