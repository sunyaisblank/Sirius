"""Seal two unexecuted original-size Alcubierre PNG actions and a prelaunch-only reseal."""
from pathlib import Path
import ast
import copy
import hashlib
import json
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
work = Path(__file__).resolve().parent
root = work.parent.parent
prior = root / 'out/forensic-kerr-calibration'
packet = root / 'out/forensic-reconstruction-b739f92'
assert Path.cwd() == root and set(p.name for p in work.iterdir()) == {'prepare.py', 'reseal.py'}


def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def write(path, value):
    assert path.is_relative_to(work) and not path.exists()
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert revision == '5733117ba8ba4430d723d62294f991db423fbf4c'
assert not subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=normal'])
base = json.loads((prior / 'preparation.json').read_text())
base_pins = json.loads((prior / 'pins.json').read_text())['files']
for name, expected in base_pins.items():
    assert identity(Path(name)) == expected, name
case = next(c for c in json.loads((packet / 'preparation.json').read_text())['cases'] if c['name'] == 'alcubierre-png')
assert case['format'] == 'png' and case['width'] == 1920 and case['height'] == 1080
assert case['requested_spp'] == case['effective_spp'] == 64 and case['rays_per_renderer'] == 132_710_400
for role in ('legacy', 'ported'):
    config = json.loads(Path(case['configs'][role]).read_text())
    output = Path(config['render']['outputPath'])
    assert not output.exists() and output.parent.is_dir() and output.suffix == '.png'
    assert case['future_render_argv_review_only'][role][0] == base['producers'][role]['path']
for name in ('run_owned.py', 'owned_guard.py', 'run_action.py', 'read_result.py'):
    shutil.copy2(prior / name, work / name)
reader = (work / 'read_result.py').read_text()
replacements = {
    'reconstruction Kerr result': 'reconstruction Alcubierre result',
    "CASE = 'legacy-kerr-original-size'": "CASE = 'legacy-alcubierre-original-size'",
    "role = preparation['role']": "role = actions[CASE]['role']",
    "{'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}": "{'seconds': 1800, 'rss_kib': 4194304, 'reserve_kib': 2097152}",
    "preparation['case']['name'] == 'kerr-png'": "preparation['case']['name'] == 'alcubierre-png' and preparation['case']['format'] == 'png'",
    "config['render']['width'] == config['render']['height'] == 512 and config['render']['samplesPerPixel'] == 32": "config['render']['width'] == 1920 and config['render']['height'] == 1080 and config['render']['samplesPerPixel'] == 64",
    "config['metric']['name'] == 'Kerr' and config['metric']['mass'] == 1 and config['metric']['spin'] == 0.9 and\n            config['observer'] == {'distance': 30.0, 'inclination': 80.0, 'azimuth': 0.0, 'fov': 60.0}":
        "config['metric']['name'] == 'Alcubierre' and config['metric']['mass'] == 1 and config['metric']['spin'] == 0 and\n            config['metric']['warpVelocity'] == 0.5 and config['metric']['bubbleRadius'] == 1 and config['metric']['bubbleSigma'] == 0.5 and\n            config['observer'] == {'distance': 20.0, 'inclination': 80.0, 'azimuth': 0.0, 'fov': 90.0}",
    r'Resolution:\s+512 x 512': r'Resolution:\s+1920 x 1080',
    r'Samples:\s+32(?: spp)?': r'Samples:\s+64(?: spp)?',
    r'Tiles:\s+64 \(spiral order\)': r'Tiles:\s+510 \(spiral order\)',
    r'Metric:\s+Kerr': r'Metric:\s+Alcubierre',
    'tiles and all(total == 64 and 0 < n <= 64 for n, total in tiles)': "(tiles or owner['stop_reason'] is not None) and all(total == 510 and 0 < n <= 510 for n, total in tiles)",
    "    result['runtime_join'] = {": "    completed_tiles = max((n for n, _ in tiles), default=0)\n    # 480 full64x64 tiles and30 final-row64x56 tiles; counter omits tile identities.\n    ray_bounds = [completed_tiles * 262144 - min(completed_tiles, 30) * 32768,\n                  completed_tiles * 262144 - max(0, completed_tiles - 480) * 32768]\n    result['runtime_join'] = {",
    "'requested_spp': 32, 'effective_spp_from_sealed_source': 25": "'requested_spp': 64, 'effective_spp_from_sealed_source': 64",
    "'completed_tile_ray_count': max(n for n, _ in tiles) * 4096 * 25,": "'completed_tile_ray_count': ray_bounds[0] if ray_bounds[0] == ray_bounds[1] else None,\n        'completed_tile_ray_bounds': ray_bounds, 'required_full_ray_count': 132710400,",
    "Declared original CPU semantics; temperature options remain unwired and selected azimuth0 remains fixed": "Declared original CPU Alcubierre lensing/background with thin disk disabled; selected azimuth0 remains fixed",
    "output = Path(preparation['output'])": "output = Path(action['output'])",
    'list(range(1, 65))': 'list(range(1, 511))',
    "comparator = ROOT / 'bin/linux-gcc/forensic-reconstruction-b739f92/comparator/image_compare'": "comparator = Path(preparation['comparator']['path'])",
    "item['path'] == str(comparator)": "item['path'] == preparation['comparator']['copied_from']",
    "'compare', 'png', '512', '512'": "'compare', 'png', '1920', '1080'",
    "'width': 512, 'height': 512": "'width': 1920, 'height': 1080",
    "'compared_pixels': 262144, 'compared_samples': 1048576": "'compared_pixels': 2073600, 'compared_samples': 8294400",
}
for before, after in replacements.items():
    assert reader.count(before) == 1, before
    reader = reader.replace(before, after)
ast.parse(reader, filename='read_result.py')
(work / 'read_result.py').write_text(reader)
runner = (work / 'run_action.py').read_text()
point = "runpy.run_path(str(work / 'run_owned.py'), run_name='__main__')"
assert runner.count(point) == 1
runner = runner.replace(point, """if sys.argv[1] == 'ported-alcubierre-original-size':
    prior_owner = json.loads((work / 'legacy-alcubierre-original-size/owner.json').read_text())
    assert prior_owner['status'] == 'terminal' and prior_owner['remaining_owned_processes'] == []
    assert not prior_owner['cleanup_errors'] and prior_owner['source_unchanged'] and prior_owner['whole_inputs_unchanged']
    assert prior_owner['observed_births'] and all(item['absent'] for item in prior_owner['observed_births'])
""" + point)
ast.parse(runner, filename='run_action.py')
(work / 'run_action.py').write_text(runner)
original_comparator = root / 'bin/linux-gcc/forensic-reconstruction-b739f92/comparator/image_compare'
receipt = json.loads((packet / 'build-comparator/receipt.json').read_text())
expected = next(item for item in receipt['products'] if item['path'] == str(original_comparator))
assert identity(original_comparator) == {k: expected[k] for k in ('bytes', 'sha256')}
shutil.copy2(original_comparator, work / 'image_compare')
assert identity(work / 'image_compare') == identity(original_comparator)
manifest = copy.deepcopy(base)
manifest.pop('role')
manifest.pop('output')
limits = {'seconds': 1800, 'rss_kib': 4194304, 'reserve_kib': 2097152}
manifest.update(case=case, limits=limits, staged_runner_revision=revision,
                outputs={role: json.loads(Path(case['configs'][role]).read_text())['render']['outputPath'] for role in ('legacy', 'ported')},
                comparator={'path': str(work / 'image_compare'), 'copied_from': str(original_comparator), **identity(work / 'image_compare')},
                budget_scope='Two original-size Alcubierre PNGs, legacy then ported;1800s per-render operational diagnostic cap only,4GiB RSS/2GiB reserve; no measured completion/performance guarantee, resized substitute or change to Vulkan300s acceptance.',
                reader_source={'path': str(prior / 'read_result.py'), **identity(prior / 'read_result.py')},
                output_contract='Full1920x1080 noninterlaced RGBA8 PNG; all2073600pixels/8294400channels under unchanged validated comparator.',
                required_render_count=2, total_rays=265_420_800, tiles_per_render=510,
                timing_evidence='No matching timing exists for the selected Alcubierre domain; observe full-scene tile progress and classify diagnostic stops as incomplete.',
                progress_scope='510tiles:480full64x64 and30edge64x56; counter-only partial completed-ray bounds preserve unknown edge membership; full count132710400.',
                live_capture={'helper': str(root / 'out/capture-forensic-live.py'), 'role_source': 'actions.json[case].role', 'scope': 'Stable-birth actual ELF, role-specific LDD libraries, cwd/argv and relevant environment only'},
                resealable_helper_paths=[str(work / 'run_action.py'), str(work / 'read_result.py'), str(root / 'out/capture-forensic-live.py')],
                input_format_scope='Public Alcubierre script emits PPM; PNG is an explicit new-reconstruction writer variant, not recovered July tape configuration')
manifest['interpretation_limits'] = [text for text in manifest['interpretation_limits'] if 'temperature' not in text]
manifest['interpretation_limits'].append('Alcubierre CPU path disables thin disk and retains original lensing/background semantics; no changed physics')
write(work / 'preparation.json', manifest)
old_action = json.loads((prior / 'actions.json').read_text())['legacy-kerr-original-size']
prefix = old_action['argv'][:-len(base['case']['future_render_argv_review_only']['legacy'])]
actions = {}
for role in ('legacy', 'ported'):
    name = role + '-alcubierre-original-size'
    actions[name] = {'revision': revision, 'role': role, 'output': manifest['outputs'][role], 'kind': 'check', 'expected_exit': 0,
                     'limits': limits, 'scope': role + ' original-size1920x1080/64SPP Alcubierre PNG only; NEW_RECONSTRUCTION,1800s operational diagnostic cap, no completion/performance/original-tape/scientific/release claim',
                     'argv': prefix + case['future_render_argv_review_only'][role]}
inputs = sorted(set(list(base_pins) + [str(prior / 'pins.json'), str(prior / 'read_result.py'), str(root / 'out/capture-forensic-live.py'),
                    str(packet / 'build-comparator/receipt.json'), str(packet / 'comparator-controls/owner.json'), str(packet / 'comparator-controls/receipt.json')] +
                    [str(p) for p in work.iterdir() if p.is_file()] + [str(work / 'actions.json'), str(work / 'pins.json')]))
for action in actions.values():
    action['frozen_inputs'] = inputs
write(work / 'actions.json', actions)
pins = {name: identity(Path(name)) for name in inputs if name != str(work / 'pins.json')}
write(work / 'pins.json', {'files': pins, 'file_count': len(pins), 'file_records_sha256': hashlib.sha256(json.dumps(pins, sort_keys=True, separators=(',', ':')).encode()).hexdigest()})
print(json.dumps({'prepared_not_executed': True, 'staged_runner_revision': revision, 'pins': len(pins),
                  'reader': identity(work / 'read_result.py'), 'preparation': identity(work / 'preparation.json'),
                  'pins_file': identity(work / 'pins.json'), 'comparator_copy': identity(work / 'image_compare'),
                  'commands': [['/usr/bin/python3.12', str(work / 'run_action.py'), name] for name in actions]}, indent=2))
