"""Gate and apply the exact three-file AP patch with the established authority."""
from pathlib import Path
from authority import source, replace
from component_gate import component, contract_inputs, doc

component()
contract_inputs(doc(Path(__file__).resolve().parent / 'provider/source-contract.json'))
source('out/native-cpu-v4-provider/apply_provider.py', 'f398377c1bf91e297342143763b4be4d11bf017b3ec3e87b1a7cc366a1b1e0b8')
body = source('out/native-cpu-normalization-tags/apply_provider.py', '21b6f4f7f86f54b65407907ce63c6e096529ae89e3e084941b1cbdb8cc7cfda2')
body = replace(body,
    'AJ caches the exact private-register normalization predicate with complete write invalidation and unchanged arithmetic/fallback',
    'AP extends only Endpoint pure/scalar admission for runtime header word11241; the unchanged shader handles acceptance/refusal and all outputs; six AO objects including merged math remain byte-identical')
exec(compile(body, str(Path(__file__)), 'exec'))
