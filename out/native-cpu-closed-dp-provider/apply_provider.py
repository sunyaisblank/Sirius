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
    'Exact AV three-file path/identity patch only. Headers, cs_exec_fn, positive job/full FP scope, first7 objects, all public branches and AU T/E/shared9 remain exact. Only admitted DP v4 calls closed DP/local cc92 math9. Compiled-interface evidence only. Original genuine-job C8/all32 is the next numerical gate; no ABBA/frame until remaining-stage closure is justified. No qualification transfer or adoption.')
exec(compile(body, str(Path(__file__)), 'exec'))
component()
