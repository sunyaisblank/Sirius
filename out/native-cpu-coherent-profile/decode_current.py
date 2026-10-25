"""Decode the existing terminal capture with explicit checkout/proof roots."""
from pathlib import Path
import hashlib
R=Path.cwd();base=R/'out/native-cpu-runtime-observe/decode_all_threads.py'
assert hashlib.sha256(base.read_bytes()).hexdigest()=='185082fa17ecdfbe5d01fb626cfc9e77b7dd440c3569d5ccb22739d75d263d97'
source=base.read_text()
changes=[("ROOT = Path(__file__).resolve().parents[2]", "ROOT = Path.cwd()"),
 ("PROOF = ROOT / 'out/native-cpu-transport-worker-occupancy/provider/linked-native-proof.json'", "PROOF = ROOT / 'out/native-cpu-local-identities/provider/linked-native-proof.json'")]
for a,b in changes:
 assert source.count(a)==1;source=source.replace(a,b,1)
exec(compile(source,str(base),'exec'))
