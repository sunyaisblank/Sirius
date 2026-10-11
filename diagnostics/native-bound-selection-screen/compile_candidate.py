"""Bounded existing toolchain screen of an isolated native result-selection candidate."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import ast, hashlib, json, os, re, signal, subprocess, time
ROOT=Path.cwd()
WORK=Path(__file__).resolve().parent
BASELINE='2053ef79c90752d37b0075b8137cfb9b571b6820'
PRIOR=ROOT/'attestations/native-vulkan/24a8e14/native-rhs-layer-clock-execution/files'
TOOLS=ROOT/'bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools'
def seal(p):
 d=p.read_bytes();return {'bytes':len(d),'sha256':hashlib.sha256(d).hexdigest()}
def read(p):return json.loads(p.read_text())
def git(*args):return subprocess.check_output(['git',*args],text=True,timeout=30).strip()
assert WORK==ROOT/'out/native-bound-selection'
assert git('rev-parse','HEAD')==BASELINE==git('rev-parse','@{u}')
assert not git('status','--porcelain')
assert read(WORK/'source-review.json')['accepted']
assert read(WORK/'compiler-controller-review.json')['accepted']
prior=read(PRIOR/'compiler.json');assert prior['accepted']
paths=[Path(__file__),WORK/'preparation.json',WORK/'prototype.patch',WORK/'source-review.json',WORK/'compiler-controller-review.json',WORK/'phase/retained_transport.slang',ROOT/'scripts/build-retained-kernels.py',PRIOR/'compiler.json']
paths+=sorted((WORK/'kernels').iterdir())
for p in (Path('/opt/slang/bin/slangc'),Path('/opt/slang/lib/libslang-compiler.so.0.2026.12.0.1'),TOOLS/'spirv-dis',TOOLS/'spirv-as',TOOLS/'spirv-val'):
 assert seal(p)==prior['input_seals'][str(p)];paths.append(p)
prep=read(WORK/'preparation.json')
for n,expected in prep['input_seals'].items():assert seal(Path(n))==expected;paths.append(Path(n))
assert seal(WORK/'kernels/retained_pair.slang')==prep['candidate_pair']
assert seal(WORK/'phase/retained_transport.slang')==prior['input_seals'][str(ROOT/'out/native-rhs-layer-clock-admission/phase/retained_transport.slang')]
record={'accepted':False,'baseline_revision':BASELINE,'input_seals':{str(p):seal(p) for p in paths},'seconds_per_command':60,'rss_bytes':536870912,'runs':[],'errors':[],'scope':'Native-default Transport result-mask compiler/validator screen; original diagnostic observer preserved. No native device/pipeline/dispatch/numerical/benefit/qualification verdict.'}
signals=[]
for sig in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):signal.signal(sig,lambda n,f:signals.append(n))
def save():
 record['signals']=list(signals);(WORK/'compiler.json').write_text(json.dumps(record,indent=2)+'\n')
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
 args=list(prior['runs'][0]['arguments']);args[1]=str(WORK/'phase/retained_transport.slang');args[args.index('-I')+1]=str(WORK/'kernels');args.insert(2,'-DSIRIUS_RETAINED_BOUND_MASK32=1')
 raw=WORK/'phase/retained_transport.compiler.spv';assembly=WORK/'phase/retained_transport.spvasm';final=WORK/'phase/retained_transport.spv';assert not any(p.exists() for p in (raw,assembly,final));args[-1]=str(raw)
 run('transport-slang',args);run('transport-disassemble',[TOOLS/'spirv-dis',raw,'-o',assembly])
 builder=ROOT/'scripts/build-retained-kernels.py';tree=ast.parse(builder.read_text());node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='native_controls');env={'re':re};exec(compile(ast.Module(body=[node],type_ignores=[]),str(builder),'exec'),env);text=env['native_controls'](assembly.read_text());assembly.write_text(text)
 run('transport-assemble',[TOOLS/'spirv-as','--target-env','spv1.5',assembly,'-o',final]);run('transport-validate',[TOOLS/'spirv-val','--target-env','vulkan1.2',final])
 assert text.count('OpReadClockKHR')==2 and text.count('OpControlBarrier')==9
 assert 'OpCapability Int64' not in text and 'OpTypeInt 64' not in text
 for fn,operation in (('RPUpAdd','OpFAdd'),('RPUpMultiply','OpFMul')):
  body=re.search(r'%'+fn+r' = OpFunction\b.*?OpFunctionEnd',text,re.S)[0]
  assert not any(op in body for op in ('OpBranchConditional','OpSelectionMerge','OpLoopMerge'))
  calls=re.findall(r'OpFunctionCall %\S+ (%\S+)',body)
  assert calls.count('%OrderedAdd32' if fn=='RPUpAdd' else '%OrderedMultiply32')==1
 record['product']={'spv':seal(final),'assembly':seal(assembly),'clock_instruction_count':2,'barriers':9,'bound_helpers_conditional_branches':0}
 assert {str(p):seal(p) for p in paths}==record['input_seals'];assert git('rev-parse','HEAD')==BASELINE and not git('status','--porcelain');assert not signals;record['accepted']=True
except BaseException as e:record['errors'].append(repr(e))
finally:save()
print(json.dumps({k:record[k] for k in ('accepted','runs','errors')}));assert record['accepted']
