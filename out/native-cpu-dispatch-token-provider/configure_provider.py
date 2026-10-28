"""Gate the original cache reconfiguration; change no options or toolchain."""
from pathlib import Path
from authority import source, replace
from component_gate import component, contract_inputs, doc

source('out/native-cpu-inline-private-provider/configure_provider.py', '9bc88314b7de6a5873c45f7980f5ed8627564eb9bb9f20fa554b8263a4f0d7e2')

source('out/native-cpu-packed-math-provider/configure_provider.py', 'b047846ec5bb08ff4ea56fe4515eb47a43c801a6de4a79ddfb2e5b3c786b398f')

source('out/native-cpu-v4-raw-provider/configure_provider.py', '727d62b94d976b306e3a370109abedead3b3942140f38d10b0d64c100c01cd12')

component()
contract_inputs(doc(Path(__file__).resolve().parent / 'provider/source-contract.json'))
source('out/native-cpu-v4-provider/configure_provider.py', '456f82bec93bccc22d78cdb941c75c1a427d5295bb579b93fa0884f7a2dd95fe')
source('out/native-cpu-endpoint-header-provider/configure_provider.py', '162a50c3e0587c60c9a925484ad4759fcdb36b9c8e0287e0a50a8a32b8ab83e1')
body = source('out/native-cpu-slang-bridge/configure_provider.py', '887a9d6f484f9791aea3ba297899b3275b84ceca726fdfd6dc83d296fef30c85')
body = replace(body, 'three sealed native objects', 'eight sealed native objects')
body = replace(body, 'link3 already-sealed native objects', 'link8 native rows, changing baseline scalar3 and the adapter in seven sealed eighth-object inputs')
exec(compile(body, str(Path(__file__)), 'exec'))
