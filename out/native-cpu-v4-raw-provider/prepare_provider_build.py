"""Gate the same dry-run inventory and original bounded build ownership."""
from pathlib import Path
from authority import source, replace
from component_gate import component, contract_inputs, doc, terminal

component()
P = Path(__file__).resolve().parent / 'provider'
contract_inputs(doc(P / 'source-contract.json'))
terminal(P / 'configure/owner.json')
source('out/native-cpu-v4-provider/prepare_provider_build.py', 'ccfd67ef343bf6c24fc141ae16b6601438c34c52063469cbef766ba59bd5718f')
source('out/native-cpu-endpoint-header-provider/prepare_provider_build.py', '5ebd23db67537967c89d562e22c43e89eb586a27411df385d1d35acbda34d4c2')
body = source('out/native-cpu-slang-bridge/prepare_provider_build.py', '5a91184fb45804da505f76552448fdb617a8eff6220909c5351bb031eee531d1')
body = replace(body, 'linking3 sealed same-source native CPU objects',
    'linking8 native rows with seven sealed inputs in the eighth and the first7 AP native objects unchanged')
exec(compile(body, str(Path(__file__)), 'exec'))
