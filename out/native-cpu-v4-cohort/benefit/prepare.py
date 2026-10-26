from pathlib import Path
import hashlib
basis=Path.cwd()/'out/native-cpu-local-identities-cohort/benefit/prepare.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='2384f8db318cc3d5623cf65066ff609a4eb5a0485f87f2704fab60441a44000c'
source=basis.read_text()
import json
acceptance=Path(__file__).resolve().parent.parent/"control-8-acceptance.json"
assert json.loads(acceptance.read_text())["optional_math"]["all_selected_optional_math_observed"]
assert source.count('out/native-cpu-local-identities/provider')==2
source=source.replace('out/native-cpu-local-identities/provider','out/native-cpu-v4-provider/provider')
assert source.count('AN local integer identity-store candidate with unchanged AM math')==1
source=source.replace('AN local integer identity-store candidate with unchanged AM math','AO separately guarded optional v4 math with exact v3 fallback and unchanged AN raw3')
exec(compile(source,str(basis),"exec"))
