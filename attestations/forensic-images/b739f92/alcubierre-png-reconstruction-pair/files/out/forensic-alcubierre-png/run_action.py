"""Root-only one-shot entry point for the sealed reconstruction preparation actions."""
from pathlib import Path
import hashlib
import json
import os
import runpy
import sys

sys.dont_write_bytecode = True
work = Path(__file__).resolve().parent
assert not sys.flags.optimize, 'ownership checks require ordinary Python'
assert len(sys.argv) == 2 and sys.argv[1] in json.loads((work / 'actions.json').read_text())
assert Path.cwd() == work.parent.parent
pins = json.loads((work / 'pins.json').read_text())
tools = json.loads((work / 'preparation.json').read_text())['tools']
for name, expected in pins['files'].items():
    path = Path(name)
    assert path.is_file(), name
    assert not path.is_symlink() or (name in tools and str(path.resolve(strict=True)) == tools[name]['resolved_path']), name
    data = path.read_bytes()
    assert len(data) == expected['bytes'] and hashlib.sha256(data).hexdigest() == expected['sha256'], name
assert not any(os.environ.get(name) for name in ('GIT_DIR', 'GIT_WORK_TREE', 'GIT_INDEX_FILE', 'GIT_OBJECT_DIRECTORY', 'GIT_ALTERNATE_OBJECT_DIRECTORIES'))
if sys.argv[1] == 'ported-alcubierre-original-size':
    prior_owner = json.loads((work / 'legacy-alcubierre-original-size/owner.json').read_text())
    assert prior_owner['status'] == 'terminal' and prior_owner['remaining_owned_processes'] == []
    assert not prior_owner['cleanup_errors'] and prior_owner['source_unchanged'] and prior_owner['whole_inputs_unchanged']
    assert prior_owner['observed_births'] and all(item['absent'] for item in prior_owner['observed_births'])
runpy.run_path(str(work / 'run_owned.py'), run_name='__main__')
