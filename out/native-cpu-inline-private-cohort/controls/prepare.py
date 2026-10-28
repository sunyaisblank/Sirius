"""Restore original C8→ABBA22/Decimal→ordered seven-control preparation."""
from pathlib import Path
import importlib.util
from control_authority import ao,load,authority_files

ao('prepare.py')
gate_path=Path(__file__).resolve().parents[1]/'benefit/correctness.py'
spec=importlib.util.spec_from_file_location('AS_current_C8_gate',gate_path)
gate=importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)
c8_inputs=gate.check()
source=load('prepare.py')
anchor='assert all((root / name).is_file() for name in inputs)'
assert source.count(anchor)==1
source=source.replace(anchor,
    'inputs |= {str(path.relative_to(root)) for path in authority_files()}\n'
    'inputs |= {str(path.relative_to(root)) for path in c8_inputs}\n'+anchor,1)
exec(compile(source,'out/native-cpu-local-identities-cohort/controls/prepare.py','exec'))
