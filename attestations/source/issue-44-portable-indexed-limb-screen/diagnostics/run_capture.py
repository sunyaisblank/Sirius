"""Bounded parent invokes one unchanged short-input guard with passive observers."""
from pathlib import Path
import hashlib
import json
import os
import resource
import subprocess
import sys

sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
resource.setrlimit(resource.RLIMIT_FSIZE,(64*1024*1024,64*1024*1024))
root=Path.cwd()
work=Path(__file__).resolve().parent
assert len(sys.argv)==2 and sys.argv[1] in ('baseline','candidate')
producer=sys.argv[1]
folder=work/('capture-'+producer)
assert folder.exists()
stage=root/'bin/linux-gcc/portable-indexed-limb-accumulation'
module=(root/'bin/linux-gcc/src/sirius/backend/retained/retained_endpoint_portable_normal_sum.spv'
        if producer=='baseline' else stage/'retained_endpoint_portable_normal_sum.spv')
observer=root/'attestations/software-vulkan/b5d2abf/optimized-ir-stage-code-binding/optimizer_observer.so'
assert hashlib.sha256(observer.read_bytes()).hexdigest()=='a1f8808cd8a3a793c59eeed10d78515a08d5aa09b1080ac47281c1a78c04fd68'
env=os.environ.copy()
env['SIRIUS_LIMB_JIT_CAPTURE']=str(folder)
script=folder/'capture.gdb'
assert not script.exists()
script.write_text('set pagination off\nset confirm off\nset breakpoint pending on\n'
                  'set debuginfod enabled off\nset print thread-events off\nset disable-randomization off\n'
                  'set environment LD_PRELOAD '+str(observer)+'\n'
                  'set environment SIRIUS_OPT_IR_FOLDER '+str(folder)+'\n'
                  'source '+str(work/'capture_short_input_gdb.py')+'\nrun\n')
record={'producer':producer,'module':str(module.relative_to(root)),
        'module_sha256':hashlib.sha256(module.read_bytes()).hexdigest(),
        'executable_sha256':hashlib.sha256((stage/'short_input_jit_probe_v2').read_bytes()).hexdigest(),
        'observer_sha256':hashlib.sha256(observer.read_bytes()).hexdigest(),
        'child_environment_additions':{'SIRIUS_LIMB_JIT_CAPTURE':str(folder),
                                      'LD_PRELOAD':str(observer),'SIRIUS_OPT_IR_FOLDER':str(folder)},
        'source_guard':'ni=1, rows=1 => ni<329; return before arithmetic and outputwrites',
        'scope':'Instrumented compiler/code-footprint observation only, no numerical or speed qualification.'}
(folder/'capture-inputs.json').write_text(json.dumps(record,indent=2)+'\n')
subprocess.run(['/usr/bin/gdb','-q','-batch','-nx','-x',str(script),'--args',
                str(stage/'short_input_jit_probe_v2'),str(module)],env=env,check=True)
state=json.loads((folder/'gdb-capture.json').read_text())
assert state['inferior'] and state['exit_code']==0 and state['error'] is None
assert any(item.get('file') for item in state['events'])
exports=folder/'optimizer-exports.tsv'
assert exports.exists() and exports.read_text().strip(), 'fresh optimizer return must be observed; cache-only dispatch is insufficient'
assert list(folder.glob('optimized-*.ll'))
(folder/'capture-verification.json').write_text(json.dumps({'passed':True,'producer':producer,
    'fresh_optimizer_exports_observed':True,'jit_registered_objects':sum('file' in item for item in state['events']),
    'scope':'Fresh compiler observation with complete unchanged outputcanary. Machine/object/IR/source joins still require separate verification; no numerical/speed/qualification claim.'},indent=2)+'\n')
