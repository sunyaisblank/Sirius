"""One ordered remaining7; stop on a failed/skipped original control."""
from pathlib import Path
import importlib.util
import json
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
runtime = work / 'runtime'
accepted = json.loads((work.parent / 'benefit/readers/benefit.json').read_text())
decimal = json.loads((work.parent / 'benefit/readers/decimal-agreement.json').read_text())
assert accepted['audit_passed'] and accepted['all22_numeric_pass']
assert decimal['agrees_on_all22'] and decimal['all22_pass']
assert not any((runtime / ('control-' + str(i))).exists() for i in range(1, 8))
spec = importlib.util.spec_from_file_location('original_controls_parser', work / 'readers/native_evidence.py')
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)
for i in range(1, 8):
    label = 'control-' + str(i)
    subprocess.run(['/usr/bin/python3', '-B', str(runtime / 'run_owned.py'), label], check=True)
    original.xml_properties(runtime / label / 'gtest.xml', label)
subprocess.run(['/usr/bin/python3', '-B', str(work / 'readers/read_controls.py')], check=True)
