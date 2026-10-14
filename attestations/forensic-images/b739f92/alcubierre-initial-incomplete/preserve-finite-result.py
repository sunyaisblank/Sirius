"""Preserve a terminal finite packet with an isolated Git index; never change HEAD."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

root = Path.cwd()
kind = sys.argv[1]
parent = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert parent == 'f5d649bf5067d1d68dc9215462f37830122f0540'
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()

def pin(path):
    raw = path.read_bytes()
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}

def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')

runner_source = parent
if kind in ('p4', 'kerr', 'volume'):
    label = {'p4': 'page-thorne', 'kerr': 'kerr-authority', 'volume': 'disk-volume'}[kind]
    case = {'p4': 'original-three-page-thorne-controls', 'kerr': 'original-four-kerr-authority-controls',
            'volume': 'original-four-disk-volume-controls'}[kind]
    work = root / ('out/science-p4-' + label)
    bundle = root / ('attestations/scientific-core/4bcd987/p4-' + label + '-current-source')
    tag = 'evidence/issue-53-p4-' + label + '-current-source'
    result = json.loads((work / 'result.json').read_text())
    assert result['passed'] and not any(result[k] for k in ('failures', 'errors', 'skips', 'disabled'))
    owner = json.loads((work / case / 'owner.json').read_text())
    assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert all(row['absent'] for row in owner['observed_births'])
    mapping = [(work, 'files/out/science-p4-' + label),
               (root / 'bin/linux-gcc/tests/oracle/sirius_oracle_tests', 'producer/sirius_oracle_tests'),
               (root / 'bin/linux-gcc/tests/oracle/sirius_oracle_tests[1]_tests.cmake', 'producer/registrations.cmake')]
    if kind in ('kerr', 'volume'):
        mapping += [(root / 'bin/linux-gcc/tests/core/sirius_core_tests', 'producer/sirius_core_tests'),
                    (root / 'bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake', 'producer/core-registrations.cmake')]
    if kind == 'volume':
        mapping += [(root / 'bin/linux-gcc/src/sirius/app/sirius_render_tests', 'producer/sirius_render_tests'),
                    (root / 'bin/linux-gcc/tests/render/sirius_render_tests[1]_tests.cmake', 'producer/render-registrations.cmake')]
    summary = result
    message = 'Preserve original current-source P4 ' + label + ' witnesses (#53)\n'
elif kind == 'native-failure':
    work = root / 'out/macos-moltenvk-integration'
    bundle = root / 'attestations/native-vulkan/4bcd987/macos-provider-discovery-f5d649b'
    tag = 'evidence/issue-48-provider-discovery-f5d649b'
    archive = work / 'native-38556777741.zip'
    assert pin(archive) == {'bytes': 75835745, 'sha256': '2b0dbcd768b646a9d258c1942c797e3a642a96ece71a4f2c2d51e2da00acf030'}
    import zipfile
    with zipfile.ZipFile(archive) as zipped:
        assert zipped.testzip() is None
        prefix = 'out/macos-moltenvk-integration/observation/'
        observation = json.loads(zipped.read(prefix + 'integration.json'))
        assert observation['runner_source']['runner_revision'] == parent
        assert observation['provider_source_inputs_unchanged'] and not observation['cases']
        assert not any(name.startswith(prefix + 'provider-build/') for name in zipped.namelist())
        owner = json.loads(zipped.read(prefix + 'toolchain/owner.json'))
        assert owner['status'] == 'terminal' and owner['child_exit'] == 1 and owner['stop_reason'] is None
        assert owner['whole_inputs_unchanged'] and owner['containment_complete'] and owner['observation_complete']
        assert not owner['cleanup_errors'] and all(row['absent'] for row in owner['observed_births'])
        assert b'CMake module directory is ambiguous' in zipped.read(prefix + 'toolchain/stderr.log')
    mapping = [(archive, 'artifact.zip'), (work / 'run-38556777741-failure.log', 'workflow-failure.log')]
    summary = {'schema': 'sirius-failed-native-provider-discovery-v1', 'run': 38556777741,
               'job': 115827784818, 'artifact': 11890208325, 'archive': pin(archive),
               'runner_revision': parent, 'provider_source_commit': 'e2422dd92017879a452377a1122f1838e6610b3c',
               'failure': 'CMake module-directory layout guessed incorrectly; query its authoritative CMAKE_ROOT instead.',
               'observation': observation, 'toolchain_owner': owner, 'provider_build_started': False,
               'qualification_claimed': False, 'raw_scope': 'Complete immutable native artifact ZIP plus source snapshots and workflow failure log. Remote-only system input pins are observations.'}
    message = 'Preserve failed native CMake discovery and complete raw artifact (#48)\n'
elif kind in ('native-escape', 'native-side-effects'):
    work = root / 'out/macos-moltenvk-integration'
    if kind == 'native-escape':
        label, runner_source, run, job, artifact = ('ninja-escape-a7a0c9d', 'a7a0c9df594da6a7944ba97c00823643a7de7075', 38557811215, 115831201718, 11890771513)
        archive_pin = {'bytes': 85241436, 'sha256': '7cd11604400b28bbe3c39e55cf1f4d86ec8f1cd64dc2aa67bcb79c0f6ca74c6d'}
    else:
        label, runner_source, run, job, artifact = ('build-side-effects-bfefab7', 'bfefab73d571b3cbeafb9f0f657924169f067239', 38559236199, 115835910513, 11891043528)
        archive_pin = {'bytes': 108476597, 'sha256': '56d46acefdbffb7e84fdcf6835eb48608d286e04d1038f99c10a47d92d5fcb6d'}
    bundle = root / ('attestations/native-vulkan/4bcd987/macos-provider-' + label)
    tag = 'evidence/issue-48-provider-' + label
    archive = work / ('native-' + str(run) + '.zip')
    assert pin(archive) == archive_pin
    import zipfile
    with zipfile.ZipFile(archive) as zipped:
        assert zipped.testzip() is None
        prefix = 'out/macos-moltenvk-integration/observation/'
        observation = json.loads(zipped.read(prefix + 'integration.json'))
        assert observation['runner_source']['runner_revision'] == runner_source
        assert not observation['cases']
        owner = json.loads(zipped.read(prefix + 'provider-build/owner.json'))
        assert owner['status'] == 'terminal'
        if kind == 'native-escape':
            assert observation['provider_source_inputs_unchanged']
            assert owner['child_exit'] == -15 and owner['stop_reason'] == 'observed descendant escaped launch group'
            assert not owner['containment_complete'] and not owner['observation_complete']
        else:
            assert not observation['provider_source_inputs_unchanged']
            assert owner['child_exit'] == 120 and owner['stop_reason'] is None
            assert owner['containment_complete'] and owner['observation_complete']
            assert len(zipped.read(prefix + 'provider-build/stderr.log')) == 64 * 1024**2
        assert owner['whole_inputs_unchanged'] and not owner['cleanup_errors']
        assert not owner['remaining_observed'] and all(row['absent'] for row in owner['observed_births'])
        for phase in ('toolchain', 'stock-inventory', 'stock-uuid'):
            prior = json.loads(zipped.read(prefix + phase + '/owner.json'))
            assert prior['child_exit'] == 0 and prior['containment_complete'] and prior['observation_complete']
            assert prior['whole_inputs_unchanged'] and not prior['cleanup_errors']
    mapping = [(work / ('run-' + str(run) + '-failure.log'), 'workflow-failure.log')]
    summary = {'schema': 'sirius-rejected-native-build-v1', 'run': run,
               'job': job, 'artifact': artifact, 'archive': pin(archive),
               'runner_revision': runner_source, 'provider_source_commit': 'e2422dd92017879a452377a1122f1838e6610b3c',
               'failure': ('Ninja compiler ABI command creates its own process group; original strict owner rejects it.' if kind == 'native-escape' else
                           'Five child-Python bytecode additions violate frozen tree; stderr exactly64MiB ends mid-dyld-line and child exits120. Truncated output does not establish a compiler semantic failure.'),
               'observation': observation, 'provider_build_owner': owner, 'qualification_claimed': False,
               'original_controls_started': False,
               'raw_scope': 'Complete native artifact ZIP split into consecutive 32MiB chunks; concatenation has the recorded archive bytes/hash. Contains complete provider source tar and all remote observations.'}
    message = 'Preserve rejected native provider build and complete raw artifact (#48)\n'
elif kind == 'forensic-incomplete':
    work = root / 'out/forensic-alcubierre-png'
    case = work / 'legacy-alcubierre-original-size'
    bundle = root / 'attestations/forensic-images/b739f92/alcubierre-initial-incomplete'
    tag = 'evidence/issue-100-alcubierre-initial-incomplete'
    result = json.loads((case / 'result.json').read_text())
    owner = json.loads((case / 'owner.json').read_text())
    assert result['classification'] == 'INCOMPLETE_DIAGNOSTIC_LIMIT' and not result['complete_original_size_render']
    assert owner['status'] == 'terminal' and owner['stop_reason'] == 'diagnostic_timeout'
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert all(row['absent'] for row in owner['observed_births'])
    assert not (work / 'ported-alcubierre-original-size').exists()
    mapping = [(work, 'files/out/forensic-alcubierre-png'),
               (root / 'bin/linux-gcc/forensic-reconstruction-b739f92/scaffold-retry/renderers/bin/Sirius', 'producer/Sirius'),
               (root / 'bin/linux-gcc/forensic-reconstruction-b739f92/scaffold-retry/renderers/bin/sirius', 'producer/sirius')]
    summary = {'schema': 'sirius-incomplete-full-size-reconstruction-v1', 'evidence_kind': 'NEW_RECONSTRUCTION',
               'runner_revision': parent, 'source_export_revision': 'b739f920c2ae67847ed7961e011d2350a839e410',
               'source_and_build_evidence_reused': 'evidence/issue-100-kerr-exr-reconstruction-pair at38c6a264f7ac18d9b75f798337eaa2027e4264ac',
               'result': result, 'owner': owner, 'paired_comparison_performed': False,
               'old_ported1800_action_launched': False, 'historical_tape_identity_claimed': False,
               'qualification_claimed': False,
               'unused_pair_helper': 'Old unused helper has stale EXR512x512 argv; no comparator was executed. The fresh completion packet must use exact PNG1920x1080 before sealing.'}
    message = 'Preserve incomplete original-size Alcubierre reconstruction without a pixel verdict (#100)\n'
else:
    raise SystemExit(kind)

assert not bundle.exists()
bundle.mkdir(parents=True)
if kind in ('native-escape', 'native-side-effects'):
    parts = bundle / 'artifact.zip.parts'
    parts.mkdir()
    hashes = []
    with archive.open('rb') as stream:
        number = 0
        while block := stream.read(32 * 1024**2):
            path = parts / f'part-{number:03d}'
            path.write_bytes(block)
            hashes.append({'path': str(path.relative_to(bundle)), **pin(path)})
            number += 1
    assert b''.join((bundle / row['path']).read_bytes() for row in hashes) == archive.read_bytes()
    summary['archive_parts'] = hashes
for source, relative in mapping:
    destination = bundle / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.is_dir():
        shutil.copytree(source, destination, symlinks=True)
    else:
        shutil.copy2(source, destination)
for name in ('scripts/macos_moltenvk_integration.py', 'scripts/macos_vulkan_identity.py',
             'scripts/diagnose-macos-retained.py', 'scripts/macos_metal_probe.py',
             'scripts/runtime_identity.py', 'scripts/governance/build_policy.py', '.github/workflows/ci.yml'):
    destination = bundle / 'runner-source' / name
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(subprocess.check_output(['git', 'show', runner_source + ':' + name]))
shutil.copy2(Path(__file__), bundle / Path(__file__).name)
write(bundle / 'result.json', summary)
write(bundle / 'manifest.json', {'qualification_claimed': False,
      'files': {str(path.relative_to(bundle)): pin(path) for path in sorted(bundle.rglob('*')) if path.is_file()}})
fd, index = tempfile.mkstemp(dir=root / 'out', prefix='preserve-finite-')
os.close(fd)
os.unlink(index)
env = dict(os.environ, GIT_INDEX_FILE=index)
try:
    subprocess.run(['git', 'read-tree', '--empty'], env=env, check=True)
    rows = []
    paths = [path for path in sorted(bundle.rglob('*')) if path.is_file()]
    for path in paths:
        assert not path.is_symlink()
        blob = subprocess.check_output(['git', 'hash-object', '-w', '--', str(path)], text=True).strip()
        rows.append(f'100644 {blob}\t{path.relative_to(root)}\0')
    subprocess.run(['git', 'update-index', '-z', '--index-info'], input=''.join(rows).encode(), env=env, check=True)
    tree = subprocess.check_output(['git', 'write-tree'], env=env, text=True).strip()
    commit = subprocess.check_output(['git', 'commit-tree', tree, '-p', runner_source], input=message.encode(), env=env).decode().strip()
    subprocess.run(['git', 'tag', tag, commit], check=True)
    for path in paths:
        raw = subprocess.check_output(['git', 'show', commit + ':' + str(path.relative_to(root))])
        assert len(raw) == path.stat().st_size and hashlib.sha256(raw).hexdigest() == pin(path)['sha256']
    print(json.dumps({'bundle': str(bundle.relative_to(root)), 'tag': tag, 'commit': commit,
                      'files': len(paths), 'manifest': pin(bundle / 'manifest.json'), 'all_git_bytes_verified': True}))
finally:
    if Path(index).exists():
        Path(index).unlink()
