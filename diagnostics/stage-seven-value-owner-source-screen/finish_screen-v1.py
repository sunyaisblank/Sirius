"""Archive the exact finished source screen and remove only its WORK.

Requires pushed useful source, specific independent clearance and exact
publication readbacks. Preserves the consumed header and all537 input records.
"""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess

ROOT = Path('/home/astra/.project/Sirius')
WORK = ROOT / 'out/stage-seven-value-owner-source-screen'
ARCHIVE = ROOT / 'attestations/review/9fbd0e8/stage7-write-mapping'
HEAD = '9fbd0e8b0d549a13476d066f2877dcd2cf9a94d9'
HEADER = ROOT / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'

def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, timeout=30)

def read(path):
    return json.loads(path.read_text())

def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

def save(path, value):
    assert not path.exists(), path
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')

def relative(path):
    return str(path.relative_to(ROOT))

def inventory():
    result = []
    assert WORK.is_dir() and not WORK.is_symlink() and WORK.resolve() == WORK
    for path in sorted(WORK.rglob('*')):
        assert not path.is_symlink(), path
        if path.is_file():
            result.append(path)
        else:
            assert path.is_dir(), path
    return result

def birth(pid):
    try:
        text = Path(f'/proc/{pid}/stat').read_text()
        fields = text[text.rfind(')') + 2:].split()
        return fields[19], int(fields[1])
    except FileNotFoundError:
        return None, None

assert Path.cwd().resolve() == ROOT
assert git('rev-parse', 'HEAD').decode().strip() == HEAD
assert not git('status', '--porcelain')
for refs in ['HEAD...@{upstream}', 'main...origin/main']:
    assert git('rev-list', '--left-right', '--count', refs).split() == [b'0', b'0']
assert git('worktree', 'list', '--porcelain').decode().count('worktree ') == 1
assert not ARCHIVE.exists()
clearance = read(WORK / 'cleanup-clearance.json')
assert clearance['accepted'] and clearance['source_revision'] == HEAD
for name, expected in clearance['approved_seals'].items():
    assert seal(WORK / name) == expected
owner = read(WORK / 'mapping/owner.json')
assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
assert owner['source']['revision'] == HEAD and owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert owner['stop_reason'] is None and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
before, after = [read(WORK / 'mapping' / f'inputs-{phase}.json') for phase in ['before', 'after']]
assert before == after and len(before) == 537
for name, expected in before.items():
    assert seal(ROOT / name) == {key: expected[key] for key in ['bytes', 'sha256']}
mapping = read(WORK / 'mapping-result.json')
assert mapping['accepted_source_mapping'] and mapping['source_revision'] == HEAD
assert seal(HEADER) == mapping['header']
assert mapping['runtime_equality_proved'] is False and mapping['performance_trial_authorized'] is False
source = read(WORK / 'source-preservation-pushed.json')
assert source['accepted'] and source['source_revision'] == HEAD and source['pushed']
assert git('ls-remote', 'origin', 'refs/tags/' + source['source_tag']).decode().split()[0] == source['source_commit']
for name, version in source['source_versions'].items():
    data = git('show', source['source_commit'] + ':' + version['git_path'])
    assert data == (WORK / name).read_bytes()
publication = read(WORK / 'publication-accepted.json')
assert publication['accepted'] and publication['source_revision'] == HEAD
for name, expected in publication['readback_seals'].items():
    assert seal(WORK / name) == expected
assert Path('/proc/sys/kernel/random/boot_id').read_text().strip() == owner['boot_id']
saved = [(item['pid'], item['saved_start_ticks']) for item in owner['observed_births']]
saved.append((owner['controller']['pid'], owner['controller']['start_ticks']))
births = []
for pid, ticks in sorted(set(saved)):
    current, _ = birth(pid)
    assert current != ticks, (pid, ticks)
    births.append({'pid': pid, 'saved_start_ticks': ticks, 'current_start_ticks': current, 'recorded_birth_absent': True})
ancestors = set()
pid = os.getpid()
while pid and pid not in ancestors:
    ancestors.add(pid)
    _, pid = birth(pid)
visible, inaccessible, matches = [], [], []
for proc in sorted(Path('/proc').iterdir()):
    if not proc.name.isdecimal() or int(proc.name) in ancestors:
        continue
    try:
        command = (proc / 'cmdline').read_bytes()
        cwd = (proc / 'cwd').resolve(strict=True)
    except (PermissionError, FileNotFoundError, ProcessLookupError) as error:
        inaccessible.append({'pid': int(proc.name), 'exception': type(error).__name__})
        continue
    visible.append(int(proc.name))
    if str(WORK).encode() in command or b'out/stage-seven-value-owner-source-screen' in command or cwd == WORK or WORK in cwd.parents:
        matches.append(int(proc.name))
assert not matches, matches
paths = inventory()
records = []
for path in paths:
    expected = seal(path)
    destination = ARCHIVE / 'files' / path.relative_to(WORK)
    destination.parent.mkdir(parents=True, exist_ok=True)
    assert not destination.exists()
    shutil.copy2(path, destination)
    assert seal(path) == expected == seal(destination)
    records.append({'origin': relative(path), 'recovery': relative(destination), **expected})
header_copy = ARCHIVE / 'consumed-header/retained_kernels.h'
header_copy.parent.mkdir(parents=True)
shutil.copy2(HEADER, header_copy)
assert seal(HEADER) == mapping['header'] == seal(header_copy)
tracked = set(git('ls-files').decode().splitlines())
recoveries = []
for name, expected in before.items():
    if name in tracked:
        data = git('show', f'{HEAD}:{name}')
        assert len(data) == expected['bytes'] and hashlib.sha256(data).hexdigest() == expected['sha256']
        recovery = f'git:{HEAD}:{name}'
    elif ROOT / name == HEADER:
        recovery = relative(header_copy)
    else:
        assert (ROOT / name).is_relative_to(WORK), name
        recovery = relative(ARCHIVE / 'files' / (ROOT / name).relative_to(WORK))
        assert seal(ROOT / recovery) == {key: expected[key] for key in ['bytes', 'sha256']}
    recoveries.append({'origin': name, 'recovery': recovery, 'bytes': expected['bytes'], 'sha256': expected['sha256']})
assert {relative(path): seal(path) for path in inventory()} == {record['origin']: {key: record[key] for key in ['bytes', 'sha256']} for record in records}
save(ARCHIVE / 'manifest.json', {'accepted': True, 'source_revision': HEAD, 'source_tag': source['source_tag'],
    'source_commit': source['source_commit'], 'records': records, 'executed_input_recoveries': recoveries,
    'header_record': {'origin': relative(HEADER), 'recovery': relative(header_copy), **mapping['header']},
    'process_check': {'saved_births': births, 'accessible_pids': visible, 'inaccessible_or_disappeared': inaccessible,
                      'excluded_ancestors': sorted(ancestors), 'WORK_matches': matches},
    'scope': 'All finished source-screen bytes preserved. Saved Linux births and accessible WORK command/cwd references only; no hidden/global/scheduled or Windows-handle claim. Source mapping only, no numeric/performance/qualification admission.'})
shutil.rmtree(WORK)
assert not WORK.exists() and seal(HEADER) == mapping['header']
assert not git('status', '--porcelain')
save(ARCHIVE / 'cleanup-result.json', {'accepted': True, 'source_revision': HEAD, 'removed': relative(WORK),
    'preserved_WORK_files': len(records), 'preserved_executed_inputs': len(recoveries), 'consumed_header_retained': True,
    'production_and_qualification_inputs_unchanged': True, 'full_goal_complete': False})
print(json.dumps({'removed': relative(WORK), 'preserved_WORK_files': len(records), 'executed_input_recoveries': len(recoveries)}))
