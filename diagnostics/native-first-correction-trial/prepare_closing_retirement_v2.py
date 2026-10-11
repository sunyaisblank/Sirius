"""Correct the unlaunched census JSON representation without changing frozen v1 files."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json
ROOT = Path.cwd(); WORK = Path(__file__).resolve().parent
assert WORK == ROOT/'out/native-division-first-correction-screen'
def seal(p):
    return {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
def put(p, value):
    assert not p.exists(), p
    p.write_text(value)
old_ps = WORK/'verify_final_retirement.ps1'; old_py = WORK/'run_final_retirement.py'
assert seal(old_ps) == {'bytes':2599,'sha256':'8631bb1590660f0937bee4940272678d0126dfc7180c763a6f6f48e6ccc06bcf'}
assert seal(old_py) == {'bytes':7542,'sha256':'b2b75419ac391117256ffe45990f75f598eba34443b357238474e685f535fe0f'}
original = old_ps.read_text(); old = ' $targets=@($births.windows_targets)\n'
new = " $targets=@($births.windows_targets|ForEach-Object {@{category=[string]$_.category;kind=[string]$_.kind;pid=[int]$_.pid;birth=[string]$_.birth}})\n"
assert original.count(old) == 1
put(WORK/'verify_final_retirement_v2.ps1', original.replace(old,new))
original = old_py.read_text()
assert original.count('verify_final_retirement.ps1') == 4
assert original.count('closing-retirement-source-review.json') == 2
new = original.replace('verify_final_retirement.ps1','verify_final_retirement_v2.ps1').replace('closing-retirement-source-review.json','closing-retirement-v2-source-review.json')
put(WORK/'run_final_retirement_v2.py', new)
put(WORK/'closing-retirement-v2-preparation.json', json.dumps({'accepted':True,
    'input_seals':{str(p.relative_to(ROOT)):seal(p) for p in (Path(__file__),old_ps,old_py)},
    'output_seals':{n:seal(WORK/n) for n in ('verify_final_retirement_v2.ps1','run_final_retirement_v2.py')},
    'supersedes_unexecuted_source_clearance':'closing-retirement-source-review.json',
    'reason':'ConvertFrom-Json produces PSCustomObjects; explicitly reconstruct mutable Hashtable targets before adding observed result key.',
    'invariants':'All51 Windows and62 Linux exact births, scoped consumers, original40s governor and ownership functions; input inventory and frozen v1 source/manifest stay unchanged.'},indent=2)+'\n')
print(json.dumps({'accepted':True,'outputs':{n:seal(WORK/n) for n in ('verify_final_retirement_v2.ps1','run_final_retirement_v2.py')}}))
