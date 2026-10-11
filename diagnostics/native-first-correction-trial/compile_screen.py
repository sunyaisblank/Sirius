"""Bounded private first-correction compiler/validator screen; no device work."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import ast, hashlib, json, os, re, signal, subprocess, time
ROOT = Path.cwd(); WORK = Path(__file__).resolve().parent
BASELINE = '2053ef79c90752d37b0075b8137cfb9b571b6820'
PRIOR = ROOT/'attestations/native-vulkan/451e2cb/native-division-phase-observer/files'
CLOCK = ROOT/'attestations/native-vulkan/2053ef7/native-shader-clock-admission/files'
ORDINARY = ROOT/'attestations/native-vulkan/3045067/ordinary-pipeline-generated-code'
COMPILER = Path('/opt/slang/bin/slangc')
LIBRARY = Path('/opt/slang/lib/libslang-compiler.so.0.2026.12.0.1')
TOOLS = ROOT/'bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):
    data=p.read_bytes();return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
def git(*args):return subprocess.check_output(['git',*args],text=True,timeout=30).strip()
assert ROOT.resolve() == Path('/home/astra/.project/Sirius')
assert WORK == ROOT/'out/native-division-first-correction-screen'
assert git('rev-parse','HEAD') == BASELINE == git('rev-parse','@{u}')
assert not git('status','--porcelain')
preparation=read(WORK/'preparation.json');baseline=read(PRIOR/'baseline-source.json')
assert preparation['accepted'] and preparation['baseline_revision'] == baseline['revision'] == BASELINE
assert preparation['source_inverse_byte_exact'] and preparation['loop_body_byte_exact']
assert set(git('ls-files').splitlines()) == set(baseline['tracked_source_seals'])
for n,v in baseline['tracked_source_seals'].items():assert seal(ROOT/n)==v,n
review_paths=set()
for name in ('source-review.json','compiler-controller-review.json'):
    review=read(WORK/name)
    assert review['accepted'] and review['baseline_revision']==BASELINE
    for path,value in review['reviewed_seals'].items():
        assert seal(Path(path))==value,path
        review_paths.add(Path(path))
assert str(Path(__file__)) in read(WORK/'compiler-controller-review.json')['reviewed_seals']
prior=read(CLOCK/'phase-compiler.json');assert prior['accepted']
for tool in (COMPILER,LIBRARY,TOOLS/'spirv-val',TOOLS/'spirv-dis',TOOLS/'spirv-as'):
    assert seal(tool)==prior['input_seals'][str(tool)]
for path,value in (preparation['input_seals']|preparation['output_seals']).items():assert seal(Path(path))==value,path
paths=set(Path(n) for n in preparation['input_seals']|preparation['output_seals'])|review_paths
paths.update(ROOT/n for n in baseline['tracked_source_seals'])
paths.update((Path(__file__),WORK/'preparation.json',WORK/'source-review.json',WORK/'compiler-controller-review.json',
              PRIOR/'compile_division_observer.py',CLOCK/'phase-compiler.json',COMPILER,LIBRARY,
              TOOLS/'spirv-val',TOOLS/'spirv-dis',TOOLS/'spirv-as'))
record={'accepted':False,'baseline_revision':BASELINE,'input_seals':{str(p):seal(p) for p in sorted(paths)},
        'seconds_per_command':60,'rss_bytes':536870912,'runs':[],'errors':[],
        'scope':'Private native-default Transport Slang/assembler/validator screen. No logical device, driver query, pipeline or dispatch; no performance/adoption/full qualification claim.'}
signals=[]
for s in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):signal.signal(s,lambda n,f:signals.append(n))
def save():
    record['signals']=list(signals)
    (WORK/'compiler.json').write_text(json.dumps(record,indent=2)+'\n')

def birth(pid):
    try:
        return Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
    except (FileNotFoundError, ProcessLookupError):
        return None

def live_group_members(group):
    members = []
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():
            continue
        try:
            fields = (p / 'stat').read_text().rsplit(')', 1)[1].split()
        except (FileNotFoundError, ProcessLookupError):
            continue
        if int(fields[2]) == group and fields[0] != 'Z':
            members.append({'pid': int(p.name), 'start_ticks': fields[19], 'state': fields[0]})
    return members

def run(name, args):
    r = {'name': name, 'arguments': list(map(str, args)), 'errors': [], 'outer_stop': False, 'peak_rss_bytes': 0}
    child = None
    ticks = None
    exit_observation = None
    started = time.monotonic()
    try:
        assert not signals
        with (WORK / (name + '.stdout')).open('xb') as out, (WORK / (name + '.stderr')).open('xb') as err:
            child = subprocess.Popen(r['arguments'], stdout=out, stderr=err, start_new_session=True)
            r.update(pid=child.pid, start_ticks=None)
            ticks = birth(child.pid)
            r['start_ticks'] = ticks
            assert ticks, 'Child birth unavailable'
            while True:
                if signals:
                    raise RuntimeError('Controller signal received')
                # Keep the exited leader unreaped until its owned group is
                # inspected and stopped. Its reserved PID prevents reuse of
                # this PGID while group cleanup is still authorized.
                status = os.waitid(os.P_PID, child.pid, os.WEXITED | os.WNOHANG | os.WNOWAIT)
                if status is not None:
                    exit_observation = {'code': status.si_code, 'status': status.si_status}
                    r['unreaped_exit_observation'] = exit_observation
                    break
                try:
                    rss = int(Path('/proc', str(child.pid), 'statm').read_text().split()[1]) * os.sysconf('SC_PAGE_SIZE')
                    r['peak_rss_bytes'] = max(r['peak_rss_bytes'], rss)
                except (FileNotFoundError, ProcessLookupError):
                    pass
                if signals or time.monotonic() - started > 60 or r['peak_rss_bytes'] > 536870912:
                    raise RuntimeError('Fixed signal/time/RSS guard triggered')
                time.sleep(.05)
            assert exit_observation == {'code': os.CLD_EXITED, 'status': 0}
    except BaseException as e:
        r['errors'].append(repr(e))
        r['outer_stop'] = True
    finally:
        if child is not None:
            leader_matches = False
            owned_unreaped = False
            try:
                # A successful waitid query also confirms this exact child
                # remains waitable by this parent, without releasing its PID.
                os.waitid(os.P_PID, child.pid, os.WEXITED | os.WNOHANG | os.WNOWAIT)
                owned_unreaped = True
                r['owned_leader_pending_reap'] = True
            except BaseException as e:
                r['errors'].append('Owned waitability observation: ' + repr(e))
            try:
                leader_matches = bool(ticks) and birth(child.pid) == ticks
                r['leader_birth_matches_before_reap'] = leader_matches
                assert leader_matches and os.getpgid(child.pid) == child.pid
            except BaseException as e:
                r['errors'].append('Unreaped leader observation: ' + repr(e))
            members = None
            try:
                members = live_group_members(child.pid)
                r['live_group_members_before_reap'] = members
            except BaseException as e:
                r['errors'].append('Pre-reap group census: ' + repr(e))
            try:
                # Close a surviving owned group regardless of leader exit.
                # Failed census also requires stop while the leader is held.
                if members is None or members:
                    assert owned_unreaped and os.getpgid(child.pid) == child.pid
                    os.killpg(child.pid, signal.SIGKILL)
                    r['forced_group_cleanup'] = True
                    if exit_observation == {'code': os.CLD_EXITED, 'status': 0}:
                        r['errors'].append('Normal rc0 required forced owned-group cleanup')
                else:
                    r['forced_group_cleanup'] = False
            except BaseException as e:
                r['errors'].append('Owned stop: ' + repr(e))
            try:
                r['returncode'] = child.wait(timeout=10)
            except BaseException as e:
                r['errors'].append('Owned reap: ' + repr(e))
            r['owned_birth_absent'] = False
            try:
                r['owned_birth_absent'] = bool(ticks) and birth(child.pid) != ticks
            except BaseException as e:
                r['errors'].append('Final birth observation: ' + repr(e))
            r['owned_group_absent'] = False
            try:
                os.killpg(child.pid, 0)
            except ProcessLookupError:
                r['owned_group_absent'] = True
            except BaseException as e:
                r['errors'].append('Group census: ' + repr(e))
        r['elapsed_seconds'] = time.monotonic() - started
        r['accepted'] = (not r['errors'] and not signals and r.get('returncode') == 0
                         and r.get('owned_birth_absent') and r.get('owned_group_absent'))
        record['runs'].append(r)
        save()
    assert r['accepted'], r

try:
    args=list(prior['reused_transport_compiler_run']['arguments'])
    expected=['-DSIRIUS_RETAINED_GENERAL_PROGRAM_WORDS=13096','-DSIRIUS_RETAINED_SCHWARZSCHILD_PROGRAM_WORDS=2871',
              '-DSIRIUS_RETAINED_SCHWARZSCHILD_LAYERS=42','-DSIRIUS_RETAINED_COEFFICIENTS=25',
              '-DSIRIUS_RETAINED_REGISTERS=459','-DSIRIUS_RETAINED_TERMS=5','-DSIRIUS_RETAINED_LANES=64',
              '-DSIRIUS_RETAINED_EXECUTION_LANES=64','-DSIRIUS_RETAINED_LAYERS=117','-DSIRIUS_RETAINED_PREFIX=0']
    assert args[2:12]==expected and args[0]==str(COMPILER)
    assert args[14:-2]==['-O0','-target','spirv','-profile','spirv_1_5','-entry','ComputeMain','-stage','compute','-denorm-mode-fp32','preserve']
    args[1]=str(WORK/'candidate/retained_transport.slang');args[args.index('-I')+1]=str(WORK/'candidate')
    products=WORK/'products';assert not products.exists();products.mkdir()
    raw=products/'transport.compiler.spv';assembly=products/'transport.spvasm';final=products/'transport.spv'
    original_assembly=products/'original-transport.spvasm'
    args[-1]=str(raw)
    run('transport-slang',args)
    run('transport-disassemble',[TOOLS/'spirv-dis',raw,'-o',assembly])
    builder=ROOT/'scripts/build-retained-kernels.py';tree=ast.parse(builder.read_text())
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='native_controls')
    env={'re':re};exec(compile(ast.Module(body=[node],type_ignores=[]),str(builder),'exec'),env)
    text=env['native_controls'](assembly.read_text());assembly.write_text(text)
    run('transport-assemble',[TOOLS/'spirv-as','--target-env','spv1.5',assembly,'-o',final])
    run('transport-validate',[TOOLS/'spirv-val','--target-env','vulkan1.2',final])
    run('original-disassemble',[TOOLS/'spirv-dis',ORDINARY/'query-payloads/transport.spv','-o',original_assembly])
    original=original_assembly.read_text()
    assert 'OpReadClock' not in text and 'ShaderClock' not in text
    assert 'OpTypeInt 64' not in text and 'OpTypeFloat 64' not in text and 'OpFma' not in text
    assert 'SignedZeroInfNanPreserve' not in text
    assert re.findall(r'OpCapability (\S+)',text)==re.findall(r'OpCapability (\S+)',original)
    assert re.findall(r'OpExecutionMode %\S+ (.*)',text)==re.findall(r'OpExecutionMode %\S+ (.*)',original)
    float_ops=lambda s:re.findall(r'= (OpF(?:Add|Sub|Mul)) ',s)
    assert float_ops(text)==float_ops(original)
    collectives=lambda s:{name:s.count(name) for name in ('OpControlBarrier','OpMemoryBarrier')}
    assert collectives(text)==collectives(original)
    def division_functions(s):
        names=dict(re.findall(r'OpName (%\S+) "([^"\n]+)"',s))
        functions={name:body for name,body in re.findall(r'(?m)^\s*(%\S+) = OpFunction .*?\n(.*?)\bOpFunctionEnd',s,re.S)}
        result=[]
        for name in names:
            if names[name]!='RTDivide':continue
            body=functions[name]
            assert len(re.findall(r'OpLoopMerge .*? DontUnroll',body))==1
            calls=[names.get(n,n) for n in re.findall(r'OpFunctionCall %\S+ (%\S+)',body)]
            result.append({'symbol':name,'helper_calls':calls,'multiply_calls':calls.count('RTMultiply'),
                           'subtract_calls':calls.count('RTSubtract'),'center_calls':calls.count('RTCenter')})
        assert len(result)==2
        return result
    before=division_functions(original);after=division_functions(text)
    assert all(f['multiply_calls']==f['subtract_calls']==2 and f['center_calls']==1 for f in before)
    assert all(f['multiply_calls']==f['subtract_calls']==3 and f['center_calls']==2 for f in after)
    record['products']={'spv':seal(final),'assembly':seal(assembly),'original_assembly':seal(original_assembly),
        'float_operation_count':len(float_ops(text)),'collectives':collectives(text),
        'original_divisions':before,'candidate_divisions':after,
        'scope':'Static helper copies and original numerical execution modes preserved; dynamic arithmetic equivalence remains source-derived and must be tested before adoption.'}
    assert {str(p):seal(p) for p in sorted(paths)}==record['input_seals']
    assert git('rev-parse','HEAD')==BASELINE==git('rev-parse','@{u}') and not git('status','--porcelain')
    assert not signals
    record['accepted']=True
except BaseException as e:record['errors'].append(repr(e))
finally:save()
print(json.dumps({k:record[k] for k in ('accepted','runs','errors')}))
assert record['accepted']
