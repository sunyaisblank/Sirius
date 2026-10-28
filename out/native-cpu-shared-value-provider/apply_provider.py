"""Apply only the three reviewed path/identity edits; preserve both AT headers."""
from pathlib import Path
from authority import source, replace
from component_gate import component, contract_inputs, doc

component()
P = Path(__file__).resolve().parent / 'provider'
contract_inputs(doc(P / 'source-contract.json'))
body = source('out/native-cpu-normalization-tags/apply_provider.py',
    '21b6f4f7f86f54b65407907ce63c6e096529ae89e3e084941b1cbdb8cc7cfda2')
body = replace(body,
    'Exact three-file provider patch and source review only. AJ caches the exact private-register normalization predicate with complete write invalidation and unchanged arithmetic/fallback; eligibility, ordered tuple/row/refusal/publication/count and FP/scheduler contracts remain. Fresh actual linked targets/full38/C8/original frame required; Z broader evidence stays predecessor evidence. No science/release claim.',
    'Exact AU three-file path/identity patch only. Both AT headers, actual positive job branch/full FP scope, first7 objects and public prefixes remain exact. New coupled scope/shared-value code has compiled-interface/static evidence only. Genuine-job C8 is the first numerical gate, followed by fresh sole ABBA22/Decimal, all8 controls and unchanged frame. No qualification transfer or adoption.')
exec(compile(body, str(Path(__file__)), 'exec'))
component()
