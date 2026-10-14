"""Prepare two original-size EXR actions from sealed reconstruction inputs; never execute them."""
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
assert Path.cwd() == root and set(p.name for p in work.iterdir()) == {'prepare.py'}


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
case = next(c for c in json.loads((packet / 'preparation.json').read_text())['cases'] if c['name'] == 'kerr-exr')
assert case['format'] == 'exr' and case['width'] == case['height'] == 512
assert case['requested_spp'] == 32 and case['effective_spp'] == 25 and case['rays_per_renderer'] == 6_553_600
for role in ('legacy', 'ported'):
    config = json.loads(Path(case['configs'][role]).read_text())
    output = Path(config['render']['outputPath'])
    assert not output.exists() and output.parent.is_dir() and output.suffix == '.exr'
    assert case['future_render_argv_review_only'][role][0] == base['producers'][role]['path']

for name in ('run_owned.py', 'owned_guard.py', 'run_action.py', 'read_result.py'):
    shutil.copy2(prior / name, work / name)
reader = (work / 'read_result.py').read_text()
replacements = {
    "CASE = 'legacy-kerr-original-size'": "CASE = 'legacy-kerr-exr-original-size'",
    "role = preparation['role']": "role = actions[CASE]['role']",
    "preparation['case']['name'] == 'kerr-png'": "preparation['case']['name'] == 'kerr-exr' and preparation['case']['format'] == 'exr'",
    "output = Path(preparation['output'])": "output = Path(action['output'])",
    "ordinary PNG output absent": "ordinary EXR output absent",
    "comparator = ROOT / 'bin/linux-gcc/forensic-reconstruction-b739f92/comparator/image_compare'": "comparator = Path(preparation['comparator']['path'])",
    "item['path'] == str(comparator)": "item['path'] == preparation['comparator']['copied_from']",
    "'compare', 'png', '512', '512'": "'compare', 'exr', '512', '512'",
    "full RGBA8 PNG": "full ZIP/HALF RGBA EXR decoded to exact binary32 words",
    "PNG self-decode refused": "EXR self-decode refused",
    "'format': 'png', 'sample_domain': 'RGBA8'": "'format': 'exr', 'sample_domain': 'decoded-RGBA-IEEE754-binary32-words-from-HALF'",
    "unexpected full PNG decode domain/result": "unexpected full EXR decode domain/result",
    "One declared-role reconstruction PNG completed and decoded": "One declared-role reconstruction EXR completed and decoded as exact binary32 words without epsilon",
}
for before, after in replacements.items():
    assert reader.count(before) == 1, before
    reader = reader.replace(before, after)
ast.parse(reader, filename='read_result.py')
(work / 'read_result.py').write_text(reader)
runner = (work / 'run_action.py').read_text()
point = "runpy.run_path(str(work / 'run_owned.py'), run_name='__main__')"
assert runner.count(point) == 1
runner = runner.replace(point, """if sys.argv[1] == 'ported-kerr-exr-original-size':
    prior_owner = json.loads((work / 'legacy-kerr-exr-original-size/owner.json').read_text())
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
manifest.update(case=case, outputs={role: json.loads(Path(case['configs'][role]).read_text())['render']['outputPath'] for role in ('legacy', 'ported')},
                comparator={'path': str(work / 'image_compare'), 'copied_from': str(original_comparator), **identity(work / 'image_compare')},
                budget_scope='Two new-reconstruction original-size Kerr EXRs, legacy then ported; each300s/4GiB/2GiB; no original July tape, resized substitute or performance qualification.',
                reader_source={'path': str(prior / 'read_result.py'), **identity(prior / 'read_result.py')},
                output_contract='Each complete512x512 linear RGBA EXR: four A/B/G/R scanline ZIP/HALF channels, full decoded RGBA IEEE754 binary32-word equality under unchanged validated TinyEXR comparator; no epsilon or PNG substitution.',
                required_render_count=2, total_rays=13_107_200,
                live_capture={'helper': str(root / 'out/capture-forensic-live.py'), 'role_source': 'actions.json[case].role',
                              'scope': 'Root captures stable-birth actual ELF, exact role-specific LDD libraries, cwd/argv and relevant environment; no full-environment/credential collection'},
                input_format_scope='Public Kerr PNG scene with explicit new-reconstruction EXR writer variant; not recovered July tape configuration')
write(work / 'preparation.json', manifest)
old_action = json.loads((prior / 'actions.json').read_text())['legacy-kerr-original-size']
prefix = old_action['argv'][:-len(base['case']['future_render_argv_review_only']['legacy'])]
actions = {}
for role in ('legacy', 'ported'):
    name = role + '-kerr-exr-original-size'
    actions[name] = {'revision': revision, 'role': role, 'output': manifest['outputs'][role], 'kind': 'check', 'expected_exit': 0,
                     'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
                     'scope': role + ' original-size512x512 Kerr EXR requested32/effective25 only; NEW_RECONSTRUCTION, no original tape/pair/performance/scientific/release claim',
                     'argv': prefix + case['future_render_argv_review_only'][role]}
inputs = sorted(set(list(base_pins) + [str(prior / 'pins.json'), str(prior / 'read_result.py'),
                                     str(packet / 'build-comparator/receipt.json'), str(packet / 'comparator-controls/owner.json'),
                                     str(packet / 'comparator-controls/receipt.json')] +
                    [str(p) for p in work.iterdir() if p.is_file()] + [str(work / 'actions.json'), str(work / 'pins.json')]))
for action in actions.values():
    action['frozen_inputs'] = inputs
write(work / 'actions.json', actions)
pins = {name: identity(Path(name)) for name in inputs if name != str(work / 'pins.json')}
write(work / 'pins.json', {'files': pins, 'file_count': len(pins),
      'file_records_sha256': hashlib.sha256(json.dumps(pins, sort_keys=True, separators=(',', ':')).encode()).hexdigest()})
print(json.dumps({'prepared_not_executed': True, 'pins': len(pins), 'reader': identity(work / 'read_result.py'),
                  'preparation': identity(work / 'preparation.json'), 'pins_file': identity(work / 'pins.json'),
                  'comparator_copy': identity(work / 'image_compare'),
                  'commands': [['/usr/bin/python3.12', str(work / 'run_action.py'), name] for name in actions]}, indent=2))
