"""Reuse the original Mesa configuration with eight native rows."""
from pathlib import Path
from authority import source, replace
from component_gate import component, contract_inputs, doc

component()
P = Path(__file__).resolve().parent / 'provider'
contract_inputs(doc(P / 'source-contract.json'))
body = source('out/native-cpu-slang-bridge/configure_provider.py',
    '887a9d6f484f9791aea3ba297899b3275b84ceca726fdfd6dc83d296fef30c85')
body = replace(body, 'three sealed native objects', 'eight pinned native objects')
body = replace(body, 'link3 already-sealed native objects',
    'link8 pinned native rows, including the twelve-input closed-DP object; both AT headers remain exact')
exec(compile(body, str(Path(__file__)), 'exec'))
