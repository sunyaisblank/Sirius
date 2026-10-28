"""Keep AO's C8 producer and seal its complete predecessor source chain."""
from rebind import ao, load, authority_files, seal_authorities

ao('prepare.py')
source = load('prepare.py')
source = source.replace('unchanged AM CPP math',
    'same original AM arithmetic with separately guarded optional v4 and original v3 fallback')
anchor = 'exec(compile(source, str(basis), "exec"))'
assert source.count(anchor) == 1
source = source.replace(anchor, 'source = seal_authorities(source)\n' + anchor, 1)
exec(compile(source, '/home/astra/.project/Sirius/out/native-cpu-local-identities-cohort/prepare.py', 'exec'))

# Describe this producer's actual object identities in future owner receipts.
# The pinned predecessor assertions, executable and frozen inputs stay exact.
import json
from pathlib import Path
_scope_path = Path.cwd() / 'out/native-cpu-closed-dp-unwind-cohort/actions.json'
_scope_actions = json.loads(_scope_path.read_text())
_scope_old = 'unchanged Transport/DopriPhase scalar objects and verified Endpoint scalar/pure admission extension at relative word11241'
_scope_new = 'unchanged AT guarded/admitted scalar3 and pure admission including Endpoint header11241; only DP unwind metadata changes; genuine positive-job DP v4 callbacks execute the same code/math/ABI with a verified full FDE; AU T/E scoped callbacks/shared9 math remain exact, with compiled-interface-only component scope and current C8/all32 as first genuine-job numerical gate'
assert _scope_actions['control-8']['scope'].count(_scope_old) == 1
_scope_actions['control-8']['scope'] = _scope_actions['control-8']['scope'].replace(_scope_old, _scope_new, 1)
_scope_path.write_text(json.dumps(_scope_actions, indent=2) + '\n')
