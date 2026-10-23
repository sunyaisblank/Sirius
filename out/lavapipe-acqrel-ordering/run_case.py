"""Original case/count/assertions on the exact existing matched ELF."""
from pathlib import Path
import importlib.util, json, subprocess, sys

R = Path.cwd(); W = Path(__file__).resolve().parent
RUN = W / 'runtime'; label, folder = sys.argv[1:]
action = json.loads((RUN / 'actions.json').read_text())[label]
spec = importlib.util.spec_from_file_location('actual_provider_binding', W / 'provider_binding.py')
binding = importlib.util.module_from_spec(spec); spec.loader.exec_module(binding); binding.admit()
for previous in action['prerequisites']:
    record = json.loads((RUN / previous / 'owner.json').read_text())
    assert record['passed'] and not record['remaining_owned_processes']
if label.startswith('control-') and label != 'control-8':
    benefit = json.loads((W / 'readers/benefit.json').read_text())
    decimal = json.loads((W / 'readers/decimal-agreement.json').read_text())
    assert benefit['audit_passed'] and benefit['all22_numeric_pass']
    assert decimal['agrees_on_all22'] and decimal['all22_pass']
subprocess.run([str(R / 'bin/linux-gcc/native-original-input-ownership/matched/audited_backend_tests'),
    '--gtest_filter=' + action['case'], '--gtest_color=no'], check=True)
