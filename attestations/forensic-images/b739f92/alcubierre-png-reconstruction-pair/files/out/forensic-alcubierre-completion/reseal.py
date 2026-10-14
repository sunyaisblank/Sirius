"""Root-only prelaunch reseal: runner revision and explicit diagnostic-helper pins only."""
from pathlib import Path
import ast
import hashlib
import json
import re
import subprocess
import sys

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parent.parent


def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def reseal():
    assert not sys.flags.optimize and Path.cwd() == ROOT
    preparation = json.loads((WORK / 'preparation.json').read_text())
    actions = json.loads((WORK / 'actions.json').read_text())
    pins = json.loads((WORK / 'pins.json').read_text())
    assert all(not (WORK / name).exists() and not Path(action['output']).exists() for name, action in actions.items()), 'reseal is prelaunch only'
    assert not subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=normal'], cwd=ROOT), 'runner must be committed and clean'
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    staged = preparation['staged_runner_revision']
    subprocess.run(['git', 'merge-base', '--is-ancestor', staged, revision], cwd=ROOT, check=True)
    changed = subprocess.check_output(['git', 'diff', '--name-only', staged, revision], cwd=ROOT, text=True).splitlines()
    assert all(name in ('.github/workflows/ci.yml', 'scripts/diagnose-macos-retained.py', 'scripts/governance/build_policy.py') or
               re.fullmatch(r'scripts/macos_[a-z0-9_]+\.py', name) for name in changed), ('runner change outside diagnostic scope', changed)
    helpers = set(preparation['resealable_helper_paths'])
    mutable_metadata = {str(WORK / 'preparation.json'), str(WORK / 'actions.json')}
    for name, saved in pins['files'].items():
        if name in helpers:
            ast.parse(Path(name).read_text(), filename=name)
        else:
            assert identity(Path(name)) == saved, ('immutable historical/ownership/input changed', name)
    previous = preparation['runner_revision']
    preparation['runner_revision'] = revision
    for action in actions.values():
        assert action['revision'] == previous and action['limits'] == {'seconds': 2400, 'rss_kib': 4194304, 'reserve_kib': 2097152}
        action['revision'] = revision
    for name, value in (('preparation.json', preparation), ('actions.json', actions)):
        temporary = WORK / (name + '.tmp')
        assert not temporary.exists()
        temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
        temporary.replace(WORK / name)
    refreshed = {}
    for name in helpers | mutable_metadata:
        assert name in pins['files']
        refreshed[name] = identity(Path(name))
        pins['files'][name] = refreshed[name]
    pins['file_records_sha256'] = hashlib.sha256(json.dumps(pins['files'], sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    temporary = WORK / 'pins.json.tmp'
    assert not temporary.exists()
    temporary.write_text(json.dumps(pins, indent=2, sort_keys=True) + '\n')
    temporary.replace(WORK / 'pins.json')
    return {'resealed': True, 'previous_runner_revision': previous, 'runner_revision': revision,
            'changed_runner_paths': changed, 'refreshed_diagnostic_or_metadata_pins': refreshed,
            'historical_source_config_producer_comparator_ownership_pins_unchanged': True,
            'pins': identity(WORK / 'pins.json')}


if __name__ == '__main__':
    print(json.dumps(reseal(), indent=2, sort_keys=True))
