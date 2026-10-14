"""Seal a fresh 2400s original-size Alcubierre pair; preserve the 1800s attempt."""
from pathlib import Path
import ast
import copy
import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
import sys

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parent.parent
PRIOR = ROOT / 'out/forensic-alcubierre-png'
RUNNER = 'f5d649bf5067d1d68dc9215462f37830122f0540'
LIMITS = {'seconds': 2400, 'rss_kib': 4194304, 'reserve_kib': 2097152}


def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def write(path, value):
    assert path.is_relative_to(WORK) and not path.exists()
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


assert Path.cwd() == ROOT and {p.name for p in WORK.iterdir()} == {'prepare.py'}
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == RUNNER
assert not subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=normal'])
old_preparation = json.loads((PRIOR / 'preparation.json').read_text())
old_actions = json.loads((PRIOR / 'actions.json').read_text())
old_pins = json.loads((PRIOR / 'pins.json').read_text())
assert old_preparation['runner_revision'] == RUNNER
assert old_preparation['limits'] == {'seconds': 1800, 'rss_kib': 4194304, 'reserve_kib': 2097152}
assert not (PRIOR / 'ported-alcubierre-original-size').exists()
old_case = PRIOR / 'legacy-alcubierre-original-size'
owner = json.loads((old_case / 'owner.json').read_text())
result = json.loads((old_case / 'result.json').read_text())
assert owner['status'] == 'terminal' and owner['stop_reason'] == 'diagnostic_timeout'
assert owner['child_exit'] == -15 and owner['remaining_owned_processes'] == [] and owner['cleanup_errors'] == []
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
assert result['classification'] == 'INCOMPLETE_DIAGNOSTIC_LIMIT'
assert not result['complete_original_size_render'] and result['output'] is None
assert not result['historical_tape_identity_claimed'] and not result['pair_pixel_identity_claimed']
join = result['runtime_join']
assert join['workers'] == 15 and join['tiles_reported'] == list(range(1, 478))
assert join['completed_tile_ray_bounds'] == [124059648, 125042688]
assert join['required_full_ray_count'] == 132710400
assert join['requested_spp'] == join['effective_spp_from_sealed_source'] == 64
assert owner['limits'] == old_preparation['limits']
for name, pin in old_pins['files'].items():
    assert identity(Path(name)) == pin, name
assert hashlib.sha256(json.dumps(old_pins['files'], sort_keys=True, separators=(',', ':')).encode()).hexdigest() == old_pins['file_records_sha256']

# Preserve every natural old file, including the old unlaunched ported action.
# These copies remain observations; none is used as an execution entry point.
snapshot = WORK / 'prior-1800s-attempt'
old_records = {}
for path in sorted(PRIOR.rglob('*')):
    assert not path.is_symlink(), path
    relative = path.relative_to(PRIOR)
    target = snapshot / relative
    info = path.stat()
    if path.is_dir():
        target.mkdir(parents=True, exist_ok=True)
        target.chmod(stat.S_IMODE(info.st_mode))
    else:
        assert path.is_file(), path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        pin = identity(path)
        assert identity(target) == pin
        assert target.stat().st_uid == info.st_uid and target.stat().st_gid == info.st_gid
        assert stat.S_IMODE(target.stat().st_mode) == stat.S_IMODE(info.st_mode)
        old_records[str(relative)] = {**pin, 'mode': stat.S_IMODE(info.st_mode),
                                      'uid': info.st_uid, 'gid': info.st_gid}
write(WORK / 'prior-attempt-manifest.json', {
    'source_root': str(PRIOR), 'snapshot_root': str(snapshot), 'classification': result['classification'],
    'scope': 'Whole-byte copy of the completed incomplete1800s attempt and original preparation; no overwrite/relabel or old ported execution.',
    'files': old_records})
for name in ('run_owned.py', 'owned_guard.py', 'run_action.py', 'image_compare'):
    shutil.copy2(PRIOR / name, WORK / name)
for name in ('read_result.py', 'reseal.py'):
    text = (PRIOR / name).read_text()
    assert text.count("{'seconds': 1800, 'rss_kib': 4194304, 'reserve_kib': 2097152}") == 1
    text = text.replace("{'seconds': 1800, 'rss_kib': 4194304, 'reserve_kib': 2097152}", repr(LIMITS))
    ast.parse(text, filename=name)
    (WORK / name).write_text(text)
assert identity(WORK / 'image_compare') == {'bytes': 397736, 'sha256': '3a224d55f8f74fb4fec0e2f377e0cea003983f2a53f01f15ce0579b2424373d6'}

# Only the new copy corrects this unused helper's format/dimension arguments.
# Its successful pair verdict still requires both unchanged full-size readers.
pair = (PRIOR / 'compare_pair.py').read_text()
assert pair.count("'compare','exr','512','512'") == 1
pair = pair.replace("'compare','exr','512','512'", "'compare','png','1920','1080'")
ast.parse(pair, filename='compare_pair.py')
(WORK / 'compare_pair.py').write_text(pair)

low, high = join['completed_tile_ray_bounds']
wall = owner['wall_seconds']
assessment = {
    'prior_status': result['classification'], 'prior_wall_seconds': wall,
    'completed_tiles': 477, 'total_tiles': 510, 'completed_tile_fraction': 477 / 510,
    'completed_ray_bounds': [low, high], 'full_rays': 132710400,
    'remaining_ray_bounds': [132710400 - high, 132710400 - low],
    'linear_completed_work_total_seconds_range': [wall * 132710400 / high, wall * 132710400 / low],
    'linear_tile_total_seconds': wall * 510 / 477,
    'selected_new_operational_seconds': 2400,
    'reason': '477 of510 completed tiles and at least93.48% completed ray work justify one finite2400s completion attempt with margin over linear completed-work extrapolation.',
    'limits': 'Tile costs and scheduling are nonuniform; extrapolation is not a completion bound or performance result. Ported timing is unmeasured. No resize, SPP reduction, checkpoint/resume, science assertion or Vulkan300s gate change.'}
manifest = copy.deepcopy(old_preparation)
manifest.update(runner_revision=RUNNER, staged_runner_revision=RUNNER, limits=LIMITS,
    budget_scope='Fresh original-size Alcubierre PNG pair, legacy then ported;2400s operational diagnostic cap per render, unchanged4GiB RSS/2GiB reserve. No performance/completion guarantee or Vulkan300s criterion change.',
    timing_evidence=assessment,
    prior_incomplete_attempt={'manifest': str(WORK / 'prior-attempt-manifest.json'),
        'snapshot': str(snapshot), 'original': str(PRIOR), 'old_ported_action_executed': False,
        'classification': result['classification'], 'complete_output': False},
    reader_source={'path': str(snapshot / 'read_result.py'), **identity(snapshot / 'read_result.py')},
    comparator={**old_preparation['comparator'], 'path': str(WORK / 'image_compare')},
    resealable_helper_paths=[str(WORK / 'run_action.py'), str(WORK / 'read_result.py'), str(ROOT / 'out/capture-forensic-live.py')])
assert manifest['evidence_kind'] == 'NEW_RECONSTRUCTION' and not manifest['historical_tape_identity_claimed']
assert manifest['case']['width'] == 1920 and manifest['case']['height'] == 1080
assert manifest['case']['requested_spp'] == manifest['case']['effective_spp'] == 64
assert manifest['case']['rays_per_renderer'] == 132710400
for role, path in manifest['outputs'].items():
    assert not Path(path).exists(), path
    config = json.loads(Path(manifest['case']['configs'][role]).read_text())
    assert config['render']['outputPath'] == path
    assert config['render']['width'] == 1920 and config['render']['height'] == 1080 and config['render']['samplesPerPixel'] == 64
write(WORK / 'preparation.json', manifest)
actions = copy.deepcopy(old_actions)
for action in actions.values():
    assert action['limits'] == old_preparation['limits']
    action['limits'] = LIMITS
    action['scope'] = action['scope'].replace('1800s operational diagnostic cap', '2400s operational diagnostic cap')
    action['revision'] = RUNNER
    del action['frozen_inputs']
inputs = {name for name in old_pins['files'] if not Path(name).is_relative_to(PRIOR)}
inputs.update(str(path) for path in WORK.rglob('*') if path.is_file())
inputs.update((str(WORK / 'actions.json'), str(WORK / 'pins.json')))
for action in actions.values():
    action['frozen_inputs'] = sorted(inputs)
write(WORK / 'actions.json', actions)
pins = {name: identity(Path(name)) for name in sorted(inputs) if name != str(WORK / 'pins.json')}
write(WORK / 'pins.json', {'files': pins, 'file_count': len(pins),
    'file_records_sha256': hashlib.sha256(json.dumps(pins, sort_keys=True, separators=(',', ':')).encode()).hexdigest()})
for name, record in old_records.items():
    assert identity(PRIOR / name) == {k: record[k] for k in ('bytes', 'sha256')}
for path in WORK.glob('*.py'):
    ast.parse(path.read_text(), filename=str(path))
print(json.dumps({'prepared_not_executed': True, 'prior_assessment': assessment,
    'snapshot_files': len(old_records), 'sealed_pins': len(pins),
    'owner_guard_byte_equal': all((WORK / name).read_bytes() == (PRIOR / name).read_bytes()
                                for name in ('run_owned.py', 'owned_guard.py')),
    'reader': identity(WORK / 'read_result.py'), 'preparation': identity(WORK / 'preparation.json'),
    'pins': identity(WORK / 'pins.json'),
    'commands': [['python3.12', '-B', str(WORK / 'run_action.py'), name] for name in actions]}, indent=2))
