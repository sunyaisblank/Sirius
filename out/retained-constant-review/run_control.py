import datetime, hashlib, json, os, pathlib, signal, subprocess, sys, time
root=pathlib.Path.cwd(); provider,case,seconds=sys.argv[1:]; seconds=float(seconds)
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
assert not subprocess.check_output(['git','status','--porcelain'],text=True)
folder=root/'out'/('native-wide-selection-'+revision[:7])/provider/case
folder.mkdir(parents=True,exist_ok=False)
env=os.environ.copy()
if provider=='dozen':
 env['VK_DRIVER_FILES']=str(root/'bin/linux-ci/driver-probes/dzn_icd.json')
 env['LD_LIBRARY_PATH']=str(root/'bin/linux-ci/driver-probes/mesa-26.2.3/usr/lib/x86_64-linux-gnu')+':/usr/lib/wsl/lib'
else:
 assert provider=='llvmpipe'; env['VK_DRIVER_FILES']='/usr/share/vulkan/icd.d/lvp_icd.json'
identity=[sys.executable,str(root/'scripts/runtime_identity.py'),'--source-root',str(root),'--source-revision',revision,'--source-tree-clean','true','--executable',str(root/'bin/linux-gcc/src/sirius/app/sirius'),'--require-vulkan','true','--output',str(folder/'identity.json')]
with (folder/'identity.log').open('w') as stream: subprocess.run(identity,env=env,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=30)
filters={
 'camera':'RetainedComputeTest.BatchedCameraPreservesPhysicalColumnsAndRejectsInvalidRows:RetainedComputeTest.SmoothRayCameraPreservesPhysicalLensDerivatives',
 'admission':'RetainedComputeAdmission.*',
 'factored':'RetainedComputeTest.FactoredEndpointsPreserveIndependentRootsAndBoundaryAdmission',
 'wide-science':'RetainedComputeTest.Fp64ProductsPreserveIndependentScienceOrDeclineUnsupportedDevices',
 'dp':'RetainedDopriTest.*',
 'connected':'RetainedComputeTest.JointRkStagesRetainCriticalIncrementsAndEmbeddedError:RetainedComputeTest.SchwarzschildStagesPreserveIndependentFieldsAndGeneralFallback:RetainedComputeTest.ProjectedEndpointsKeepPhysicalColumnsAndRetainedContinuation:RetainedComputeTest.DenseSegmentsPreserveSmallCovariantArrivalDerivatives:RetainedComputeTest.PhysicalInitializationRetainsTheHamiltonianResidual:RetainedComputeTest.CoupledIntervalsRequireEmbeddedAndIndependentDenseAgreement:RetainedComputeTest.SharedEndpointDensePreservesPrivateIntervalsAndSerialRetry:RetainedComputeTest.SharedTracerCompletesDeviceIntervalsAndRetainsRollbackState:RetainedComputeTest.RejectedStepRowsCannotExposeOldOrPartialCandidates',
 'timestamps':'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults',
}
exe=root/'bin/linux-gcc/tests/backend/sirius_backend_tests'
cmd=[str(exe),'--gtest_filter='+filters[case],'--gtest_color=no','--gtest_output=xml:'+str(folder/'gtest.xml')]
def birth(pid):
 try: return (pathlib.Path('/proc')/str(pid)/'stat').read_text().split(') ',1)[1].split()[19]
 except FileNotFoundError: return None
utc=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
record={'source_revision':revision,'source_tree_clean_at_start':True,'provider':provider,'case':case,'start_utc':utc(),'owner_pid':os.getpid(),'owner_birth_ticks':birth(os.getpid()),'test_executable_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'command':cmd,'outer_guard_seconds':seconds,'rss_guard_bytes':4*1024**3,'outer_stop':False,'status':'running','scope':'Original finite controls; diagnostic host RSS/time guard invalidates incomplete trials. Product allocations, numerical budgets and references unchanged; no cold/frame/full qualification claim.'}
write=lambda name,obj:(folder/name).write_text(json.dumps(obj,indent=2)+'\n')
samples=[]; start=time.monotonic(); child=None; pidfd=None
try:
 with (folder/'stdout.log').open('wb') as output:
  child=subprocess.Popen(cmd,env=env,stdout=output,stderr=subprocess.STDOUT,start_new_session=True)
  record.update(test_pid=child.pid,test_birth_ticks=birth(child.pid)); pidfd=os.pidfd_open(child.pid); write('owner.json',record)
  while child.poll() is None:
   status=pathlib.Path('/proc')/str(child.pid)/'status'
   fields={}
   try:
    for line in status.read_text().splitlines():
     key,_,value=line.partition(':')
     if key in ['State','VmRSS','VmHWM','VmSwap']: fields[key]=value.strip()
   except FileNotFoundError: pass
   rss=int(fields.get('VmRSS','0 kB').split()[0])*1024
   samples.append({'elapsed_seconds':time.monotonic()-start,'rss_bytes':rss,**fields})
   reason='host_rss' if rss>record['rss_guard_bytes'] else ('elapsed' if time.monotonic()-start>seconds else None)
   if reason:
    record.update(outer_stop=True,guard_trigger=reason)
    signal.pidfd_send_signal(pidfd,signal.SIGCONT); signal.pidfd_send_signal(pidfd,signal.SIGTERM)
    try: child.wait(timeout=5)
    except subprocess.TimeoutExpired: signal.pidfd_send_signal(pidfd,signal.SIGKILL); child.wait(timeout=5)
    break
   time.sleep(1)
  code=child.wait(timeout=5)
finally:
 if child is not None and child.poll() is None:
  signal.pidfd_send_signal(pidfd,signal.SIGCONT); signal.pidfd_send_signal(pidfd,signal.SIGKILL); child.wait(timeout=5)
 if pidfd is not None: os.close(pidfd)
 record.update(status='invalidated' if record['outer_stop'] else 'completed',returncode=child.returncode if child else None,elapsed_seconds=time.monotonic()-start,end_utc=utc(),peak_sampled_rss_bytes=max((s['rss_bytes'] for s in samples),default=0),source_tree_clean_at_end=not subprocess.check_output(['git','status','--porcelain'],text=True),owned_process_absent=birth(child.pid)!=record.get('test_birth_ticks') if child else True)
 write('memory.json',samples); write('owner.json',record)
print(json.dumps({k:record[k] for k in ['source_revision','provider','case','status','returncode','elapsed_seconds','peak_sampled_rss_bytes','owned_process_absent']}),flush=True)
sys.exit(0 if record['returncode']==0 and not record['outer_stop'] else 1)
