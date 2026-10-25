"""Execute the sealed packing-adapter compile once."""
from pathlib import Path
import hashlib
import json
import subprocess

R=Path.cwd()
W=Path(__file__).resolve().parent
prepared=json.loads((W/'prepared-bridge.json').read_text())


def pin(path):
    return {'bytes':path.stat().st_size,'sha256':hashlib.file_digest(path.open('rb'),'sha256').hexdigest()}


def verify():
    assert not subprocess.check_output(['git','status','--porcelain']).strip()
    assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==prepared['revision']
    for name,expected in prepared['inputs'].items():
        assert pin(R/name)==expected,name


verify()
output=R/prepared['object']
assert not output.exists() and not (W/'bridge-build-result.json').exists()
output.parent.mkdir(parents=True,exist_ok=True)
code=subprocess.run(prepared['argv'],cwd=R).returncode
result={'passed':code==0,'argv':prepared['argv'],'exit':code,'new_CPP_compilations':1,
        'arithmetic_compilations':0,'scope':prepared['scope']}
if output.exists():result['object']={'path':prepared['object'],**pin(output)}
verify()
(W/'bridge-build-result.json').write_text(json.dumps(result,indent=2)+'\n')
raise SystemExit(code)
