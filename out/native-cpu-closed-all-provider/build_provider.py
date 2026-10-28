"""Build/stage once with the original limits; do not load the provider."""
from pathlib import Path
from authority import source, replace
from component_gate import component, contract_inputs, doc, terminal

component()
P = Path(__file__).resolve().parent / 'provider'
contract_inputs(doc(P / 'source-contract.json'))
terminal(P / 'configure/owner.json')
body = source('out/native-cpu-slang-bridge/build_provider.py',
    'c3d87afecb5ecb975c4dc49d8cf7a30f54e07cc553418d57e6f967b66b3b6325')
body = replace(body, 'native-cpu-slang-bridge/provider', 'native-cpu-closed-all-provider/provider')
exec(compile(body, str(Path(__file__)), 'exec'))
