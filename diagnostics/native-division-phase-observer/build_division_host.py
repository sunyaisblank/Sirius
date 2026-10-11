"""Owned existing-preset source build; no default all/Mandatory target."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,os,signal,subprocess,time
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
assert WORK==ROOT/'out/native-division-phase-observer'
assert ROOT.resolve()==Path('/home/astra/.project/Sirius')
def git(*args):return subprocess.check_output(['git',*args],text=True,timeout=30).strip()
def seal(p):
 b=p.read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def read(p):return json.loads(p.read_text())
candidate=read(WORK/'candidate-source.json');HEAD=candidate['source_revision']
assert candidate['baseline_revision']=='2053ef79c90752d37b0075b8137cfb9b571b6820'
assert git('rev-parse','HEAD')==HEAD==git('rev-parse','@{u}') and not git('status','--porcelain')
assert git('branch','--show-current')==candidate['development_branch']
changed={'src/sirius/backend/vulkan/vulkan_device.cpp','tests/backend/retained_compute_test.cpp','tests/backend/shader_clock_observer.h'}
assert set(candidate['changed_source_seals'])==changed
assert git('rev-parse','HEAD^')==candidate['baseline_revision']
assert set(git('diff-tree','--no-commit-id','--name-only','-r','HEAD').splitlines())==changed
review_paths=set()
for name,required in (('host-source-review.json',('src/sirius/backend/vulkan/vulkan_device.cpp','tests/backend/retained_compute_test.cpp','tests/backend/shader_clock_observer.h')),
                      ('build-controller-review.json',())):
 review=read(WORK/name);assert review['accepted'] and review['baseline_revision']==candidate['baseline_revision']
 assert all(str(ROOT/p) in review['reviewed_seals'] for p in required)
 if name=='build-controller-review.json':
  assert str(Path(__file__)) in review['reviewed_seals'] and str(WORK/'candidate-source.json') in review['reviewed_seals']
 for p,v in review['reviewed_seals'].items():assert seal(Path(p))==v;review_paths.add(Path(p))
for p,v in candidate['changed_source_seals'].items():assert seal(ROOT/p)==v

paths=[ROOT/name for name in git('ls-files').splitlines()]
paths += [Path(__file__),WORK/'host-source-review.json',WORK/'build-controller-review.json',WORK/'candidate-source.json']
paths=sorted(set(paths)|review_paths)
record={'accepted':False,'source_revision':HEAD,'input_seals':{str(p):seal(p) for p in paths},'runs':[],'errors':[],'seconds_per_command':300,'sampled_owned_group_rss_bytes':4294967296,'scope':'Existing linux-gcc preset qualification-mode source/product build only. All original source/scientific/runtime controls remain required.'}
signals=[]
for sig in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):signal.signal(sig,lambda n,f:signals.append(n))
def save():
 record['signals']=list(signals);(WORK/'source-build.json').write_text(json.dumps(record,indent=2)+'\n')
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
                    rss = 0
                    for member in live_group_members(child.pid):
                        try:
                            rss += int(Path('/proc', str(member['pid']), 'statm').read_text().split()[1]) * os.sysconf('SC_PAGE_SIZE')
                        except (FileNotFoundError, ProcessLookupError):
                            pass
                    r['peak_rss_bytes'] = max(r['peak_rss_bytes'], rss)
                except (FileNotFoundError, ProcessLookupError):
                    pass
                if signals or time.monotonic() - started > 300 or r['peak_rss_bytes'] > 4294967296:
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
 run('configure',['cmake','--preset','linux-gcc','-DSIRIUS_ALIGNMENT_MODE=qualification'])
 cache=(ROOT/'bin/linux-gcc/CMakeCache.txt').read_text()
 assert 'SIRIUS_ALIGNMENT_MODE:STRING=qualification' in cache and 'SIRIUS_MANDATORY_TESTS:BOOL=ON' in cache and 'SIRIUS_REQUIRE_VULKAN_RUNTIME:BOOL=ON' in cache
 run('build',['cmake','--build','--preset','linux-gcc','--target','SiriusAlignmentGate','SiriusSourceGovernance','sirius','sirius_backend_tests','-j2'])
 assert {str(p):seal(p) for p in paths}==record['input_seals']
 assert git('rev-parse','HEAD')==HEAD==git('rev-parse','@{u}') and not git('status','--porcelain')
 assert not signals;record['accepted']=True
except BaseException as e:record['errors'].append(repr(e))
finally:save()
print(json.dumps({k:record[k] for k in ('accepted','source_revision','runs','errors')}));assert record['accepted']
