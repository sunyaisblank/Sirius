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
    'Exact AX three-file path/identity patch only. Headers, cs_exec_fn, positive-job/full FPscope, first7 objects and all public/scalar/ISA fallback branches remain exact. Private T/E v4 select closed canonical entries; verified closed DP unwind object and local-effects math9 are reused. Eleven constituents, single9-core family. Source/compiler/static evidence gives no numerical transfer. Fresh original current C8/all32 is required before all original downstream gates.')
exec(compile(body, str(Path(__file__)), 'exec'))
component()
