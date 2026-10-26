"""Gate the original cache reconfiguration; change no options or toolchain."""
from pathlib import Path
from authority import source, replace
from component_gate import component, contract_inputs, doc

component()
contract_inputs(doc(Path(__file__).resolve().parent / 'provider/source-contract.json'))
source('out/native-cpu-v4-provider/configure_provider.py', '456f82bec93bccc22d78cdb941c75c1a427d5295bb579b93fa0884f7a2dd95fe')
body = source('out/native-cpu-slang-bridge/configure_provider.py', '887a9d6f484f9791aea3ba297899b3275b84ceca726fdfd6dc83d296fef30c85')
body = replace(body, 'three sealed native objects', 'eight sealed native objects')
body = replace(body, 'link3 already-sealed native objects', 'link8 sealed objects, replacing only two Endpoint routing objects')
exec(compile(body, str(Path(__file__)), 'exec'))
