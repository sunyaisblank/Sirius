#!/usr/bin/env python3
"""One supervised grouped diagnostic, using default Lavapipe options."""
import datetime
import ctypes
import hashlib
import json
import os
from pathlib import Path
import signal
import struct
import subprocess
import time
import xml.etree.ElementTree as ET
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
REPORT=HERE/'run-report.json'
TIMEOUT=600
RSS_GUARD_KIB=12*1024*1024
RESERVE_KIB=2*1024*1024
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(r):
    temporary=REPORT.with_suffix('.tmp');temporary.write_text(json.dumps(r,indent=2)+'\n');temporary.replace(REPORT)
def rss_tree(pid):
    pending=[pid];seen=set();total=0
    while pending:
        current=pending.pop()
        if current in seen:continue
        seen.add(current)
        try:
            lines=Path(f'/proc/{current}/status').read_text().splitlines()
            total+=sum(int(line.split()[1]) for line in lines if line.startswith('VmRSS:'))
            pending.extend(map(int,Path(f'/proc/{current}/task/{current}/children').read_text().split()))
        except (FileNotFoundError,ProcessLookupError):pass
    return total
def available():
    return next(int(line.split()[1]) for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:'))
def start_ticks(pid):
    return int(Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()[19])
def interrupted(signum, _frame):
    raise InterruptedError(f'supervisor received signal {signum}')
for signum in (signal.SIGTERM,signal.SIGINT):signal.signal(signum,interrupted)
assert not REPORT.exists(),'one grouped run only; do not relaunch'
preparation=json.loads((HERE/'preparation-report.json').read_text())
for name,value in preparation['bound_diagnostic_files'].items():assert sha(HERE/name)==value,name
for name,value in preparation['production_sources_sha256'].items():assert sha(ROOT/name)==value,name
for field in ['candidate_sources_sha256','candidate_artifacts_sha256','link_library_sha256']:
    for name,value in preparation[field].items():assert sha(ROOT/name)==value,name
assert sha(HERE.parent/'stage-artifact-report.json')==preparation['stage_report_sha256']
assert sha(HERE.parent/'baseline.json')==preparation['baseline_sha256']
assert sha(HERE.parent/'verify_artifacts.py')==preparation['offline_verifier_sha256']
driver_options={k:v for k,v in os.environ.items() if k.startswith(('GALLIVM_','LP_','MESA_'))}
assert not driver_options,driver_options
assert not os.environ.get('VK_DRIVER_FILES') and not os.environ.get('VK_ADD_DRIVER_FILES')
icd=Path('/usr/share/vulkan/icd.d/lvp_icd.json');icd_contents=json.loads(icd.read_text())
soname=icd_contents['ICD']['library_path']
loaded_driver=ctypes.CDLL(soname)
mapped_libraries={Path(line.split(maxsplit=5)[5]).resolve(strict=True)
                  for line in Path('/proc/self/maps').read_text().splitlines()
                  if len(line.split(maxsplit=5))==6 and Path(line.split(maxsplit=5)[5]).name==Path(soname).name}
assert len(mapped_libraries)==1,mapped_libraries
library=mapped_libraries.pop()
assert library.is_file()
assert available()>=RESERVE_KIB,'insufficient initial available-memory reserve'
environment=os.environ.copy();environment['VK_ICD_FILENAMES']=str(icd);environment['SIRIUS_VULKAN_DEVICE']='0'
evidence=HERE/'evidence';evidence.mkdir(exist_ok=True)
xml_path=evidence/'transport.xml'
assert not xml_path.exists() and not (evidence/'fp32-products').exists() and not (evidence/'fp64-products').exists()
command=[str(HERE/'transport_runner'),str(HERE/'modules'),str(evidence),'--gtest_color=no',f'--gtest_output=xml:{xml_path}',
         '--gtest_filter=BothProductModes/TransportPrototype.JointRkStagesRetainCriticalIncrementsAndEmbeddedError/*']
r={'scope':preparation['scope'],'status':'prepared','started_utc':now(),'command':command,
   'source_head':preparation['source_head'],'preparation_sha256':sha(HERE/'preparation-report.json'),'runner_source_sha256':sha(Path(__file__)),
   'supervisor_pid':os.getpid(),'supervisor_start_ticks':start_ticks(os.getpid()),
   'boot_id':Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
   'concurrent_activity':'Native Radeon P1 invocation completed before this launch; other host work is possible. This is a feasibility/resource observation, not isolated performance measurement.',
   'executable_sha256':sha(HERE/'transport_runner'),'bounds':{'timeout_seconds':TIMEOUT,'sampled_process_tree_rss_guard_kib':RSS_GUARD_KIB,
                                                            'available_memory_reserve_kib':RESERVE_KIB},
   'driver_options':{'GALLIVM_':'unset','LP_':'unset','MESA_':'unset'},
   'selection_environment':{'VK_ICD_FILENAMES':str(icd),'SIRIUS_VULKAN_DEVICE':'0'},
   'icd':{'path':str(icd),'sha256':sha(icd),'contents':icd_contents},
   'driver_library':{'path':str(library),'sha256':sha(library),'bytes':library.stat().st_size,
                     'resolution':'Actual ctypes.CDLL ICD SONAME mapping in supervisor /proc/self/maps before launch; runner independently reports actual resident driver path.'},
   'peak_sampled_process_tree_rss_kib':0,'minimum_sampled_available_kib':available(),
   'guard_reason':None,'terminated_by_supervisor':False,'terminal_exit_code':None}
save(r);started=time.monotonic()
with (HERE/'raw.log').open('w') as stdout,(HERE/'stderr.log').open('w') as stderr:
    process=subprocess.Popen(command,cwd=ROOT,env=environment,stdout=stdout,stderr=stderr,start_new_session=True)
    r['pid']=process.pid
    try:
        r['process_start_ticks']=start_ticks(process.pid)
        assert os.getpgid(process.pid)==process.pid and os.getsid(process.pid)==process.pid
        r['mapped_executable_sha256']=sha(Path(f'/proc/{process.pid}/exe'))
        assert r['mapped_executable_sha256']==r['executable_sha256']
        r['status']='running';save(r)
        while process.poll() is None:
            elapsed=time.monotonic()-started;rss=rss_tree(process.pid);avail=available()
            r['peak_sampled_process_tree_rss_kib']=max(r['peak_sampled_process_tree_rss_kib'],rss)
            r['minimum_sampled_available_kib']=min(r['minimum_sampled_available_kib'],avail)
            reason='timeout' if elapsed>TIMEOUT else 'process_rss_guard' if rss>RSS_GUARD_KIB else 'available_memory_reserve' if avail<RESERVE_KIB else None
            if reason:
                r['guard_reason']=reason;r['terminated_by_supervisor']=True
                os.killpg(process.pid,signal.SIGTERM)
                try:process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGKILL);process.wait(timeout=5)
                break
            time.sleep(0.1)
        r['terminal_exit_code']=process.wait(timeout=5)
    except BaseException as error:
        r['supervisor_error']=repr(error)
        if process.poll() is None:
            assert start_ticks(process.pid)==r['process_start_ticks'] and os.getpgid(process.pid)==process.pid
            r['terminated_by_supervisor']=True;os.killpg(process.pid,signal.SIGKILL)
        r['terminal_exit_code']=process.wait(timeout=5)
    finally:
        r['whole_child_wall_seconds']=time.monotonic()-started;r['finished_utc']=now();r['status']='terminal';save(r)
observations=[json.loads(line) for line in (HERE/'raw.log').read_text().splitlines() if line.startswith('{')]
r['observations']=observations
r['raw_log_sha256']=sha(HERE/'raw.log');r['stderr_sha256']=sha(HERE/'stderr.log')
r['diagnostic_files_unchanged']=all(sha(HERE/name)==value for name,value in preparation['bound_diagnostic_files'].items())
r['production_sources_unchanged']=all(sha(ROOT/name)==value for name,value in preparation['production_sources_sha256'].items())
r['candidate_sources_unchanged']=all(sha(ROOT/name)==value for name,value in preparation['candidate_sources_sha256'].items())
r['candidate_artifacts_unchanged']=all(sha(ROOT/name)==value for name,value in preparation['candidate_artifacts_sha256'].items())
r['link_libraries_unchanged']=all(sha(ROOT/name)==value for name,value in preparation['link_library_sha256'].items())
r['offline_identity_proofs_unchanged']=(sha(HERE.parent/'stage-artifact-report.json')==preparation['stage_report_sha256'] and
                                       sha(HERE.parent/'baseline.json')==preparation['baseline_sha256'] and
                                       sha(HERE.parent/'verify_artifacts.py')==preparation['offline_verifier_sha256'])
r['driver_library_unchanged']=sha(library)==r['driver_library']['sha256']
r['driver_library_after']={'path':str(library),'sha256':sha(library),'bytes':library.stat().st_size}
r['icd_unchanged']=sha(icd)==r['icd']['sha256']
r['no_live_owned_process']=not Path(f'/proc/{process.pid}').exists()
identities=[v for v in observations if v.get('device_identity')==1]
terminal=[v for v in observations if v.get('mode_terminal')==1]
r['actual_resident_driver_matches']=bool(identities) and all(Path(v['resident_driver_path']).resolve()==library.resolve() and sha(Path(v['resident_driver_path']))==r['driver_library']['sha256'] for v in identities)
r['mode_evidence']={}
for label in ['fp32','fp64']:
    directory=evidence/(label+'-products');mode={}
    originals={kind:sha(directory/f'loaded-original-{kind}.spv') for kind in ['Camera','Transport','Endpoint','Dense','Initialize','RayCamera'] if (directory/f'loaded-original-{kind}.spv').exists()}
    mode['loaded_originals_sha256']=originals
    mode['all_six_original_modules_exactly_matched']=len(originals)==6 and all(value==preparation['original_modules'][label+':'+kind]['sha256'] for kind,value in originals.items())
    output=directory/'scientific-output.bin'
    if output.exists():
        data=output.read_bytes();mode['scientific_output_sha256']=sha(output);mode['scientific_output_bytes']=len(data)
        mode['scientific_output_shape_complete']=len(data)==48132 and struct.unpack_from('<III',data)==(0x50545352,1,15)
    readback=directory/'production-readback.bin'
    if readback.exists():mode['production_readback_sha256']=sha(readback);mode['production_readback_bytes']=readback.stat().st_size
    r['mode_evidence'][label]=mode
xml_tests=[]
if xml_path.exists():
    root=ET.parse(xml_path).getroot();r['terminal_xml_sha256']=sha(xml_path)
    for test in root.iter('testcase'):
        xml_tests.append({'name':test.attrib.get('name'),'status':test.attrib.get('status'),'result':test.attrib.get('result'),
                          'seconds':test.attrib.get('time'),'failures':len(test.findall('failure')),'skipped':len(test.findall('skipped')),
                          'properties':{p.attrib['name']:p.attrib['value'] for p in test.findall('./properties/property')}})
r['terminal_tests']=xml_tests
r['pass']=(r['terminal_exit_code']==0 and not r['terminated_by_supervisor'] and len(xml_tests)==2 and
           all(v['status']=='run' and v['result']=='completed' and v['failures']==v['skipped']==0 for v in xml_tests) and
           len(identities)==len(terminal)==2 and {v['product_mode'] for v in terminal}=={'fp32','fp64'} and
           all(v['fixture_contract_pass']==1 and v['dispatch_attempts']==v['completed_dispatches']==1 and
               v['allocation_count']==12 and v['required_buffer_bytes']<=v['resident_bytes']<=8*1024*1024 and
               v['stage_submissions']==[0,1,0,0,0,0] for v in terminal) and
           all(v['all_six_original_modules_exactly_matched'] and v.get('scientific_output_shape_complete') for v in r['mode_evidence'].values()) and
           r['actual_resident_driver_matches'] and r['diagnostic_files_unchanged'] and r['production_sources_unchanged'] and
           r['candidate_sources_unchanged'] and r['candidate_artifacts_unchanged'] and r['link_libraries_unchanged'] and
           r['offline_identity_proofs_unchanged'] and r['driver_library_unchanged'] and r['icd_unchanged'] and r['no_live_owned_process'])
r['disposition']='finite_transport_PASS' if r['pass'] else 'incomplete_no_numerical_verdict' if r['terminated_by_supervisor'] or not xml_tests else 'terminal_test_failure'
save(r)
print(json.dumps({'disposition':r['disposition'],'terminal_exit_code':r['terminal_exit_code'],'wall_seconds':r['whole_child_wall_seconds'],
                  'peak_rss_kib':r['peak_sampled_process_tree_rss_kib'],'guard_reason':r['guard_reason'],'mode_evidence':r['mode_evidence']}))
raise SystemExit(0 if r['pass'] else 1)
