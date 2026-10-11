"""Compile only two native default controls and their private guarded candidates."""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import re
import struct
import subprocess
import sys
import time
sys.dont_write_bytecode = True
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
commands = []

def seal(path):
    path = Path(path)
    raw = path.read_bytes()
    return {'path':str(path),'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}

def proc(pid):
    fields = Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()
    return {'pid':pid,'birth':fields[19],'group':int(fields[2]),'session':int(fields[3])}

def save():
    (WORK/'compiler-commands.json').write_text(json.dumps(commands,indent=2)+'\n')

def observed_run(argv, **kwargs):
    assert kwargs == {'check':True}, kwargs
    index = len(commands)
    started = time.monotonic()
    record = {'argv':argv,'tool':seal(argv[0]),'status':'prepared',
              'scope':'Tool inherits encompassing owned session/RSS/reserve guard; no independent session.'}
    commands.append(record)
    outpath, errpath = WORK/f'compiler-{index:02d}.stdout', WORK/f'compiler-{index:02d}.stderr'
    child = None
    try:
        with outpath.open('wb') as out, errpath.open('wb') as err:
            child = subprocess.Popen(argv,stdout=out,stderr=err)
            record.update(process=proc(child.pid),status='running')
            assert record['process']['session']==os.getsid(0)
            save()
            code = child.wait(timeout=290)
            record.update(status='terminal',exit=code)
            if code:
                raise subprocess.CalledProcessError(code,argv)
    finally:
        if child is not None and child.poll() is None:
            record['forced_cleanup'] = True
            child.terminate()
            try:
                child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                child.kill();child.wait(timeout=5)
        record['seconds'] = time.monotonic()-started
        record['leader_reaped'] = child is not None and child.returncode is not None
        if child is not None:
            record['exit'] = child.returncode
        for key,path in (('stdout',outpath),('stderr',errpath)):
            if path.exists():record[key]=seal(path)
        save()
    return subprocess.CompletedProcess(argv,code)

def import_file(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

emission=json.loads((WORK/'emission.json').read_text())
assert emission['pass'] and emission['production_adopted'] is False
builder=import_file('private_narrow_builder',WORK/'builder.py')
programs=import_file('original_retained_program',ROOT/'src/sirius/kernels/retained_program.py')
cache=(ROOT/'bin/linux-gcc/CMakeCache.txt').read_text()
tools={name:re.search(r'^SIRIUS_'+name+r':FILEPATH=(.*)$',cache,re.M)[1]
       for name in ('SLANGC','SPIRV_AS','SPIRV_DIS','SPIRV_VAL','SPIRV_OPT')}
tracked_header=ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
inputs=[Path(__file__).resolve(),WORK/'emit.py',WORK/'builder.py',WORK/'emission.json',
        WORK/'model.py',WORK/'model-result.json',WORK/'model-independent-review.json',
        ROOT/'scripts/build-retained-kernels.py',tracked_header,ROOT/'bin/linux-gcc/CMakeCache.txt']
inputs+=sorted(p for p in (WORK/'kernels').iterdir() if p.is_file())
inputs+=sorted(p for p in (ROOT/'src/sirius/kernels').iterdir() if p.is_file())
inputs+=[Path(path) for path in tools.values()]
for kind in ('Transport','Endpoint'):
    inputs.append(ROOT/f'bin/linux-gcc/src/sirius/backend/retained/retained_{kind.lower()}.spv')
before=[seal(p) for p in inputs]
(WORK/'compiler-inputs-before.json').write_text(json.dumps(before,indent=2)+'\n')
subprocess.run=observed_run
results=[]
started=time.monotonic()
for kind in ('Transport','Endpoint'):
    stem='retained_'+kind.lower()
    program=getattr(programs,'build_'+kind.lower()+'_program')(parallel=True)
    special=getattr(programs,'build_schwarzschild_'+kind.lower()+'_program')(parallel=True)
    root_count=len(program['outputs'])
    encoded=[program['instructions'],program['registers'],root_count]+program['outputs']+program['operations']+program['layer_offsets']
    special_encoded=[special['instructions'],special['registers'],root_count]+special['outputs']+special['operations']+special['layer_offsets']
    table=re.search(r'k'+kind+r'Program\{\{(.*?)\}\};',tracked_header.read_text(),re.S)[1]
    header_words=[int(x,0) for x in re.findall(r'(0x[0-9a-fA-F]+|\d+)u?',table)]
    assert encoded+special_encoded==header_words
    frame=(len(encoded),len(special_encoded),len(special['layer_offsets'])-1)
    if kind=='Endpoint':frame+= (special['prefix_instructions'],)
    options={kind.lower()+'_specialized':frame}
    outputs={}
    for name,fma in (('baseline',False),('candidate',True)):
        destination=WORK/f'{name}-{stem}.spv'
        code=builder.compile_shader(WORK/'kernels'/f'{stem}.slang',destination,
            tools['SLANGC'],tools['SPIRV_AS'],tools['SPIRV_DIS'],tools['SPIRV_VAL'],
            program['registers'],5,len(program['layer_offsets'])-1,
            program.get('prefix_instructions',0),fp64=False,portable=False,
            optimizer=tools['SPIRV_OPT'],fma=fma,**options)
        assert destination.read_bytes()==struct.pack('<'+'I'*len(code),*code)
        observed_run([tools['SPIRV_DIS'],str(destination),'-o',str(destination.with_suffix('.spvasm'))],check=True)
        if name=='baseline':
            assert destination.read_bytes()==(ROOT/f'bin/linux-gcc/src/sirius/backend/retained/{stem}.spv').read_bytes(), 'default module does not reproduce exact current product'
        outputs[name]={'module':seal(destination),'assembly':seal(destination.with_suffix('.spvasm'))}
    results.append({'kind':kind,'registers':program['registers'],'roots':root_count,
                    'general_nodes':program['instructions'],'special_nodes':special['instructions'],
                    'general_layers':len(program['layer_offsets'])-1,
                    'special_layers':len(special['layer_offsets'])-1,'specialized_frame':frame,
                    'original_program_words_exact':True,'baseline_module_exact':True,'outputs':outputs})
after=[seal(p) for p in inputs]
(WORK/'compiler-inputs-after.json').write_text(json.dumps(after,indent=2)+'\n')
assert before==after
assert len(commands)==20 and all(c['status']=='terminal' and c['exit']==0 and c['leader_reaped'] and not c.get('forced_cleanup') for c in commands)
result={'pass':True,'production_adopted':False,'source_revision':emission['source_revision'],
        'device_workloads':0,'commands':len(commands),'whole_input_count':len(before),
        'inputs_unchanged':True,'wall_seconds':time.monotonic()-started,'stages':results,
        'scope':'Two ordinary O0 SPIR-V1.5 native default controls reproduce exact current modules; two guarded FMA private candidates assemble and validate under Vulkan1.2. No device/control/benefit verdict.'}
(WORK/'compiler-result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
