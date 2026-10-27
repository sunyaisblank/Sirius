"""Gate the existing incremental build and stage only the private AR provider."""
from pathlib import Path
from authority import source, replace
from component_gate import component, contract_inputs, doc, terminal

source('out/native-cpu-v4-raw-provider/build_provider.py', '9eb14a0239348203ab16101d158607fb152fdc1c49048ab30384bcf734366f0b')

component()
P = Path(__file__).resolve().parent / 'provider'
contract_inputs(doc(P / 'source-contract.json'))
terminal(P / 'configure/owner.json')
source('out/native-cpu-v4-provider/build_provider.py', '13ce720faa86c7ada921812a7b0d771af4eb774c5a236260a0ac43e63bc374a9')
source('out/native-cpu-endpoint-header-provider/build_provider.py', '7a4abf4ad614887e6959688d554c154ea82a8905dc4f3758d766f1067252edce')
body = source('out/native-cpu-slang-bridge/build_provider.py', 'c3d87afecb5ecb975c4dc49d8cf7a30f54e07cc553418d57e6f967b66b3b6325')
body = replace(body, 'native-cpu-slang-bridge/provider', 'native-cpu-packed-math-provider/provider')
exec(compile(body, str(Path(__file__)), 'exec'))
