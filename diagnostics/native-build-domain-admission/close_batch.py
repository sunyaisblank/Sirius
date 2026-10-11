"""Preserve this exact finite batch; remove only its finished WORK on review.

Canonical qualification inputs and current build products remain in place.
Run prepare before publication/review; run finish after pushed source/readbacks.
The process check is scoped to saved births and accessible Linux references.
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
WORK = ROOT / 'out/native-build-domain-admission'
ARCHIVE = ROOT / 'attestations/alignment/9fbd0e8/native-build-domain-admission'
CANONICAL = ROOT / 'attestations/native-build/9fbd0e8'
HEAD = '9fbd0e8b0d549a13476d066f2877dcd2cf9a94d9'
MAIN = '5e706073a3fd2f9afcf8efbcde2f094122050f6a'
CASES = ['exports', 'verify-windows', 'verify-macos', 'configure', 'build',
         'check-alignment', 'reject-incomplete-release']

def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, timeout=30)

def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

def read(path):
    return json.loads(path.read_text())

def write_new(path, value):
    assert not path.exists(), path
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')

def relative(path):
    return str(path.relative_to(ROOT))

def files(directory):
    assert directory.is_dir() and not directory.is_symlink(), directory
    result = []
    for path in sorted(directory.rglob('*')):
        assert not path.is_symlink(), path
        if path.is_file():
            result.append(path)
        else:
            assert path.is_dir(), path
    return result

def copy_exact(origin, destination):
    before = seal(origin)
    if destination.exists():
        assert seal(destination) == before, destination
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(origin, destination)
    assert seal(origin) == before == seal(destination)
    return {'origin': relative(origin), 'recovery': relative(destination), **before}

def common():
    assert Path.cwd().resolve() == ROOT and WORK.resolve() == WORK
    assert git('rev-parse', 'HEAD').decode().strip() == HEAD
    assert git('rev-parse', 'main').decode().strip() == MAIN
    assert not git('status', '--porcelain')
    assert git('rev-list', '--left-right', '--count', 'HEAD...@{upstream}').split() == [b'0', b'0']
    assert git('rev-list', '--left-right', '--count', 'main...origin/main').split() == [b'0', b'0']
    assert git('worktree', 'list', '--porcelain').decode().count('worktree ') == 1
    manifest = read(ARCHIVE / 'manifest.json')
    assert manifest['source_revision'] == HEAD
    assert len(manifest['records']) == 123 and len(manifest['canonical_bundle_records']) == 80
    for record in manifest['records'] + manifest['canonical_bundle_records']:
        expected = {key: record[key] for key in ['bytes', 'sha256']}
        assert seal(ROOT / record['origin']) == expected
        assert seal(ROOT / record['recovery']) == expected
    assert len(files(CANONICAL)) == 80
    result = read(WORK / 'root-batch-result.json')
    assert result['accepted'] and result['source_revision'] == HEAD
    assert result['admitted_domains'] == ['macos-native-build', 'windows-native-build']
    assert len(result['pending_domains']) == 6 and result['actual_ready'] is False
    assert result['actual_release_ready'] is False and result['ultimate_ideal_satisfied'] is False
    identities = {}
    for case in CASES:
        owner = read(WORK / case / 'owner.json')
        assert owner['passed'] and owner['status'] == 'terminal'
        assert owner['source']['revision'] == HEAD and owner['source_unchanged']
        assert owner['whole_inputs_unchanged'] and owner['stop_reason'] is None
        assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['child_exit'] == (1 if case == 'reject-incomplete-release' else 0)
        before, after = [read(WORK / case / f'inputs-{phase}.json') for phase in ['before', 'after']]
        assert before == after
        for name, value in before.items():
            assert name not in identities or identities[name] == value
            identities[name] = value
            assert seal(ROOT / name) == {key: value[key] for key in ['bytes', 'sha256']}
    assert len(identities) == 623
    assert not (WORK / 'rejected-release-receipt.json').exists()
    return manifest, identities

def prepare():
    manifest, identities = common()
    tracked = set(git('ls-files').decode().splitlines())
    recovery = []
    for name, value in sorted(identities.items()):
        expected = {key: value[key] for key in ['bytes', 'sha256']}
        if name in tracked:
            data = git('show', f'{HEAD}:{name}')
            assert expected == {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
            recovery.append({'origin': name, 'recovery': f'git:{HEAD}:{name}', **expected})
        elif name.startswith(relative(WORK) + '/'):
            recovery.append(copy_exact(ROOT / name, ARCHIVE / 'files' / (ROOT / name).relative_to(WORK)))
        elif name.startswith(relative(CANONICAL) + '/'):
            recovery.append(copy_exact(ROOT / name, ARCHIVE / 'native-build-bundles' / (ROOT / name).relative_to(CANONICAL)))
        else:
            raise AssertionError(('unclassified executed input', name))
    review = read(WORK / 'local-consumer-review.json')
    volatile = []
    for product in review['products'].values():
        path = ROOT / product['path']
        assert seal(path) == product['executable_seal']
        volatile.append(copy_exact(path, ARCHIVE / 'volatile-products' / path.relative_to(ROOT / 'bin/linux-gcc')))
    header = ROOT / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
    assert seal(header) == review['production_header']
    volatile.append(copy_exact(header, ARCHIVE / 'volatile-products' / header.relative_to(ROOT / 'bin/linux-gcc')))
    for name in ['generated/sirius/alignment_receipt.json', 'src/sirius/app/resources/model/alignment_receipt.json', 'CMakeCache.txt']:
        path = ROOT / 'bin/linux-gcc' / name
        volatile.append(copy_exact(path, ARCHIVE / 'volatile-products' / name))
    assert seal(ROOT / 'bin/linux-gcc/generated/sirius/alignment_receipt.json') == seal(WORK / 'generated-alignment-receipt.json')
    assert seal(ROOT / 'bin/linux-gcc/src/sirius/app/resources/model/alignment_receipt.json') == seal(WORK / 'staged-app-alignment-receipt.json')
    assert seal(ROOT / 'bin/linux-gcc/CMakeCache.txt') == seal(WORK / 'final-CMakeCache.txt')
    write_new(WORK / 'executed-input-recovery.json', {'accepted': True, 'source_revision': HEAD,
              'input_records': recovery, 'volatile_records': volatile, 'canonical_bundles_retained': True,
              'initial_manifest_seal': seal(ARCHIVE / 'manifest.json'),
              'scope': 'All623 recorded executed inputs recoverable from immutable Git or exact local bytes; six volatile product/receipt/cache files archived. Raw evidence remains local.'})
    print(json.dumps({'prepared': True, 'executed_inputs': len(recovery), 'volatile_files': len(volatile)}))

def birth(pid):
    try:
        text = Path(f'/proc/{pid}/stat').read_text()
        fields = text[text.rfind(')') + 2:].split()
        return fields[19], int(fields[1])
    except FileNotFoundError:
        return None, None

def scoped_process_check():
    boot = Path('/proc/sys/kernel/random/boot_id').read_text().strip()
    saved = {}
    for case in CASES:
        owner = read(WORK / case / 'owner.json')
        assert owner['boot_id'] == boot
        for item in owner['observed_births']:
            saved[(item['pid'], item['saved_start_ticks'])] = 'observed child'
        controller = owner['controller']
        saved[(controller['pid'], controller['start_ticks'])] = 'controller'
    observation = read(WORK / 'readiness-observation.json')
    saved[(observation['pid'], observation['start_ticks'])] = 'readiness CLI'
    checked = []
    for (pid, ticks), role in sorted(saved.items()):
        current, _ = birth(pid)
        assert current != ticks, (pid, ticks, role)
        checked.append({'pid': pid, 'saved_start_ticks': ticks, 'current_start_ticks': current,
                        'recorded_birth_absent': True, 'role': role})
    ancestors = set()
    pid = os.getpid()
    while pid and pid not in ancestors:
        ancestors.add(pid)
        _, pid = birth(pid)
    visible, hidden, matches = [], [], []
    for proc in sorted(Path('/proc').iterdir()):
        if not proc.name.isdecimal() or int(proc.name) in ancestors:
            continue
        try:
            command = (proc / 'cmdline').read_bytes()
            cwd = (proc / 'cwd').resolve(strict=True)
        except (PermissionError, FileNotFoundError, ProcessLookupError) as error:
            hidden.append({'pid': int(proc.name), 'exception': type(error).__name__})
            continue
        visible.append(int(proc.name))
        if str(WORK).encode() in command or b'out/native-build-domain-admission' in command or cwd == WORK or WORK in cwd.parents:
            matches.append({'pid': int(proc.name), 'cwd': str(cwd)})
    assert not matches, matches
    return {'boot_id': boot, 'saved_births': checked, 'visible_pids': visible,
            'unreadable_or_disappeared': hidden, 'excluded_ancestors': sorted(ancestors),
            'WORK_reference_matches': matches,
            'scope': 'Fresh saved Linux births and accessible command/cwd references only; no hidden/global/scheduled consumer or Windows handle claim.'}

def finish():
    common()
    recovery = read(WORK / 'executed-input-recovery.json')
    assert recovery['accepted']
    for record in recovery['input_records'] + recovery['volatile_records']:
        if not record['recovery'].startswith('git:'):
            assert seal(ROOT / record['recovery']) == {key: record[key] for key in ['bytes', 'sha256']}
    for filename in ['source-preservation-pushed.json', 'closing-source-preservation.json',
                     'publication-accepted.json', 'cleanup-clearance.json']:
        receipt = read(WORK / filename)
        assert receipt['accepted'] and receipt['source_revision'] == HEAD
    clearance = read(WORK / 'cleanup-clearance.json')
    assert clearance['approved_initial_manifest'] == seal(ARCHIVE / 'manifest.json')
    assert clearance['approved_helper'] == seal(WORK / 'close_batch.py')
    assert clearance['approved_recovery'] == seal(WORK / 'executed-input-recovery.json')
    for filename in ['source-preservation-pushed.json', 'closing-source-preservation.json']:
        receipt = read(WORK / filename)
        assert git('ls-remote', 'origin', 'refs/tags/' + receipt['source_tag']).decode().split()[0] == receipt['source_commit']
        for name, version in receipt['source_versions'].items():
            data = git('show', receipt['source_commit'] + ':' + version['git_path'])
            assert seal(WORK / name) == {key: version[key] for key in ['bytes', 'sha256']}
            assert len(data) == version['bytes'] and hashlib.sha256(data).hexdigest() == version['sha256']
    references = scoped_process_check()
    records = [copy_exact(path, ARCHIVE / 'files' / path.relative_to(WORK)) for path in files(WORK)]
    expected = {record['origin']: {key: record[key] for key in ['bytes', 'sha256']} for record in records}
    assert {relative(path): seal(path) for path in files(WORK)} == expected
    # The only deletion target is this canonical, symlink-free, fully preserved WORK.
    write_new(ARCHIVE / 'closing-before-removal.json', {'accepted': True, 'source_revision': HEAD,
              'records': records, 'process_check': references, 'canonical_bundles_retained': True,
              'scope': 'All finite-batch WORK bytes preserved; final actual deletion result follows separately. Full goal remains unfinished.'})
    shutil.rmtree(WORK)
    assert not WORK.exists() and len(files(CANONICAL)) == 80
    write_new(ARCHIVE / 'cleanup-result.json', {'accepted': True, 'source_revision': HEAD,
              'removed': relative(WORK), 'removed_files': len(records), 'canonical_bundle_files_retained': 80,
              'production_products_retained': True, 'working_tree_clean': not bool(git('status', '--porcelain')),
              'process_scope': references['scope'], 'full_qualification_complete': False})
    print(json.dumps({'removed': relative(WORK), 'preserved_files': len(records), 'canonical_files_retained': 80}))

if __name__ == '__main__':
    assert len(sys.argv) == 2 and sys.argv[1] in ['prepare', 'finish']
    {'prepare': prepare, 'finish': finish}[sys.argv[1]]()
