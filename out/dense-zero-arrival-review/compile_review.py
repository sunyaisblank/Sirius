from pathlib import Path
import hashlib
import importlib.util
import json
import os
import signal
import subprocess
import sys
import time

sys.dont_write_bytecode = True
r = Path.cwd()
w = Path(__file__).resolve().parent
def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def dump(path, value): path.write_text(json.dumps(value, indent=2) + '\n')
def birth(pid): return Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
def members(group):
    found = []
    for p in Path('/proc').iterdir():
        if not p.name.isdecimal(): continue
        try:
            fields = (p/'stat').read_text().rsplit(')', 1)[1].split()
            if int(fields[2]) == group and int(fields[3]) == group:
                found.append({'pid':int(p.name), 'start_ticks':fields[19]})
        except (FileNotFoundError, ProcessLookupError): pass
    return found
def resources(group):
    rss = 0
    for entry in members(group):
        try:
            status = Path('/proc', str(entry['pid']), 'status').read_text()
            rss += next((int(x.split()[1]) for x in status.splitlines() if x.startswith('VmRSS:')), 0)
        except (FileNotFoundError, ProcessLookupError): pass
    return rss

generator = module('retained_program', r/'src/sirius/kernels/retained_program.py')
builder = module('retained_builder', r/'scripts/build-retained-kernels.py')
program = generator.build_dense_program(parallel=True)
words = [program['instructions'], program['registers'], 40] + program['outputs'] + program['operations'] + program['layer_offsets']
toolroot = r/'bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools'
toolpaths = ['/opt/slang/bin/slangc', str(toolroot/'spirv-as'), str(toolroot/'spirv-dis'), str(toolroot/'spirv-val'), '/usr/bin/spirv-opt']
variants = [('native',False,False), ('native-wide',True,False), ('portable',False,True), ('portable-wide',True,True)]
if len(sys.argv) == 3 and sys.argv[1] == '--variant':
    name, wide, portable = next(x for x in variants if x[0] == sys.argv[2])
    builder.compile_shader(r/'src/sirius/kernels/retained_dense.slang', w/(name+'.spv'),
        *toolpaths[:4], program['registers'], 5, len(program['layer_offsets'])-1,
        fp64=wide, portable=portable, optimizer=toolpaths[4])
    sys.exit(0)

assert len(sys.argv) == 1
identity = [f'static const uint kDenseCanonicalProgramWords = {len(words)}u;',
            f'static const uint kDenseCanonicalProgram[{len(words)}] = {{']
identity.extend(','.join(str(v)+'u' for v in words[i:i+16])+',' for i in range(0,len(words),16))
identity.append('};')
(w/'retained_dense_program_identity.slang').write_text('\n'.join(identity)+'\n')
import re
header = r/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
array = re.search(r'kDenseProgram\{\{(.*?)\}\};', header.read_text(), re.S)
assert array and list(map(int,re.findall(r'(\d+)u',array[1]))) == words
source = [Path(__file__).resolve(), r/'scripts/build-retained-kernels.py',
          *sorted((r/'src/sirius/kernels').glob('retained*')), header,
          w/'retained_dense_program_identity.slang', *map(Path,toolpaths)]
seals = {str(p):{'bytes':p.stat().st_size,'sha256':digest(p)} for p in source}
records = []
for name, wide, portable in variants:
    result = {'variant':name, 'source_head':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
              'source_diff_sha256':hashlib.sha256(subprocess.check_output(['git','diff','--binary'])).hexdigest(),
              'whole_input_seals':seals, 'program_words':len(words), 'program_exact_to_accepted_host_array':True,
              'compiler_timeout_seconds':90, 'maximum_group_rss_mib':4096, 'sample_cadence_seconds':.25,
              'controller_pid':os.getpid(), 'controller_start_ticks':birth(os.getpid()),
              'peak_sampled_rss_kib':0, 'stop_reason':None}
    args = [sys.executable,'-B',str(Path(__file__).resolve()),'--variant',name]
    result['arguments'] = args
    started = time.monotonic()
    births = {}
    with (w/(name+'.stdout')).open('wb') as out, (w/(name+'.stderr')).open('wb') as err:
        child = subprocess.Popen(args,stdout=out,stderr=err,start_new_session=True)
        result.update(child_pid=child.pid,child_start_ticks=birth(child.pid))
        try:
            while child.poll() is None:
                births.update((x['pid'],x['start_ticks']) for x in members(child.pid))
                result['peak_sampled_rss_kib'] = max(result['peak_sampled_rss_kib'],resources(child.pid))
                if time.monotonic()-started > 90: result['stop_reason']='compiler_time_bound'
                if result['peak_sampled_rss_kib'] > 4096*1024: result['stop_reason']='compiler_rss_bound'
                if result['stop_reason']: break
                time.sleep(.25)
        finally:
            if child.poll() is None:
                os.killpg(child.pid,signal.SIGTERM)
                try: child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid,signal.SIGKILL); child.wait(timeout=10)
    result.update(returncode=child.returncode,elapsed_seconds=time.monotonic()-started,
                  known_births=[{'pid':pid,'start_ticks':ticks} for pid,ticks in sorted(births.items())],
                  owned_group_remaining=members(child.pid))
    result['unchanged_inputs'] = all(p.stat().st_size == seals[str(p)]['bytes'] and digest(p) == seals[str(p)]['sha256'] for p in source)
    binary = w/(name+'.spv')
    if binary.exists(): result['compiled_module']={'bytes':binary.stat().st_size,'sha256':digest(binary)}
    result['pass_'] = child.returncode == 0 and result['stop_reason'] is None and not result['owned_group_remaining'] and result['unchanged_inputs']
    dump(w/(name+'-compile.json'),result)
    records.append(result)
    print(json.dumps({k:result[k] for k in ('variant','returncode','elapsed_seconds','peak_sampled_rss_kib','stop_reason','pass_')}),flush=True)
    if not result['pass_']: break
dump(w/'compile-results.json',{'pass_':len(records)==4 and all(x['pass_'] for x in records),'variants':records})
sys.exit(0 if len(records)==4 and all(x['pass_'] for x in records) else 1)
