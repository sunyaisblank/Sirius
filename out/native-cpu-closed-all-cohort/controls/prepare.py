"""Require current AX C8 and sole original ABBA22/Decimal before preparation."""
from pathlib import Path
import importlib.util
from control_authority import ao, load, authority_files

ao('prepare.py')
gate_path = Path(__file__).resolve().parents[1] / 'benefit/correctness.py'
spec = importlib.util.spec_from_file_location('AX_current_C8_gate', gate_path)
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)
c8_inputs = gate.check() + [gate_path]
source = load('prepare.py')
anchor = 'assert all((root / name).is_file() for name in inputs)'
assert source.count(anchor) == 1
source = source.replace(anchor,
    'inputs |= {str(path.relative_to(root)) for path in authority_files()}\n'
    'inputs |= {str(path.relative_to(root)) for path in c8_inputs}\n' + anchor, 1)
scope = 'Same genuine H host/AXprovider, originalcases/inputs/allocations/cache/governors/limits/assertions.'
assert source.count(scope) == 1
source = source.replace(scope, scope +
    ' Current AX closed3/local9 static composition and genuine-job C8/all32 are bound;'
    ' C8 directly observed only Endpoint/DopriPhase SIMD and Transport12 scalar.'
    ' No Transport SIMD execution or prior numerical/performance acceptance transfers.', 1)
exec(compile(source, 'out/native-cpu-local-identities-cohort/controls/prepare.py', 'exec'))
