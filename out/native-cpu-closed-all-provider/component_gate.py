"""Expose only the owned AX11-part compile-interface closure."""
from pathlib import Path
import hashlib
import importlib.util
import json
R=Path.cwd();C=R/'out/native-cpu-closed-all-component';path=C/'component.py'
expected=json.loads((C/'build-plan.json').read_text())['inputs'][str(path.relative_to(R))]
assert path.stat().st_size==expected['bytes'] and hashlib.sha256(path.read_bytes()).hexdigest()==expected['sha256']
spec=importlib.util.spec_from_file_location('AX_owned_closed_all_component',path)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
component,doc,pin,checked=module.component,module.doc,module.pin,module.checked
terminal,contract_inputs=module.terminal,module.contract_inputs
HEADERS,CORES,LOCAL_CORES,DP_CALLS=module.HEADERS,module.CORES,module.LOCAL_CORES,module.DP_CALLS
