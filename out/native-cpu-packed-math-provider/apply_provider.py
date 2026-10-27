"""Gate and apply the exact three-file AR patch with the established authority."""
from pathlib import Path
from authority import source, replace
from component_gate import component, contract_inputs, doc

source('out/native-cpu-v4-raw-provider/apply_provider.py', 'c088a76275dc252eb0f21ce06183fcac17c3c1890d0241d3b35ec891d266f99f')

component()
contract_inputs(doc(Path(__file__).resolve().parent / 'provider/source-contract.json'))
source('out/native-cpu-v4-provider/apply_provider.py', 'f398377c1bf91e297342143763b4be4d11bf017b3ec3e87b1a7cc366a1b1e0b8')
source('out/native-cpu-endpoint-header-provider/apply_provider.py', '20b6d6880c430dd0452d89f146b5c8179f4d88fd1506c295f033f53a9730b733')
body = source('out/native-cpu-normalization-tags/apply_provider.py', '21b6f4f7f86f54b65407907ce63c6e096529ae89e3e084941b1cbdb8cc7cfda2')
body = replace(body,
    'AJ caches the exact private-register normalization predicate with complete write invalidation and unchanged arithmetic/fallback',
    'AR composes the sealed packed math-boundary raw3 with separately qualified zero-residual v4 arithmetic; first7 AP native objects and original raw/scalar fallback remain byte-identical')
exec(compile(body, str(Path(__file__)), 'exec'))
