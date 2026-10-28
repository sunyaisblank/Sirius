"""Expose only the owned merged AW repair and unchanged AV source authorities."""
from pathlib import Path
import hashlib
import importlib.util
import json

R = Path.cwd()
C = R / 'out/native-cpu-closed-dp-unwind-component'
path = C / 'component.py'
expected = json.loads((C / 'build-plan.json').read_text())['inputs'][str(path.relative_to(R))]
assert path.stat().st_size == expected['bytes'] and hashlib.sha256(path.read_bytes()).hexdigest() == expected['sha256']
spec = importlib.util.spec_from_file_location('AW_owned_unwind_component', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
component, doc, pin, checked = module.component, module.doc, module.pin, module.checked
terminal, contract_inputs = module.terminal, module.contract_inputs
HEADERS, CORES, LOCAL_CORES, DP_CALLS = module.HEADERS, module.CORES, module.LOCAL_CORES, module.DP_CALLS
