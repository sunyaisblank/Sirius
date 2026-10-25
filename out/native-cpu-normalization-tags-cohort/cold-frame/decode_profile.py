"""Reuse the existing terminal-only decoder with actual AJ code associations."""
from pathlib import Path
import hashlib
R=Path.cwd();W=Path(__file__).resolve().parent
p=R/'out/native-cpu-runtime-observe/decode_all_threads.py';source=p.read_text()
a="PROOF = ROOT / 'out/native-cpu-transport-worker-occupancy/provider/linked-native-proof.json'";b="PROOF = ROOT / 'out/native-cpu-normalization-tags/provider/linked-native-proof.json'"
assert source.count(a)==1
exec(compile(source.replace(a,b),str(p),'exec'))
