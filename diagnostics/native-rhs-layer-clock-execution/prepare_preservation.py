"""Adapt the accepted complete-source/raw recovery workflow to this finite batch."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import ast, hashlib, json
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
OLD = ROOT / 'attestations/native-vulkan/8f1118d/native-shader-clock-execution/files/preserve_batch.py'
source = OLD.read_text()
text = source.replace('native-shader-clock-execution', 'native-rhs-layer-clock-execution')
text = text.replace('out/native-rhs-layer-clock-execution', 'out/native-rhs-layer-clock-admission')
text = text.replace('8f1118dd56d10cdd0f932613427b68409a1f77f6', '24a8e145d1da5dca75a4fa26ec8e639020eecf37')
text = text.replace('8f1118d', '24a8e14')
text = text.replace("WORK / 'workflow-final.json'", "WORK / 'ci.json'")
text = text.replace("len(r['targets']) == 6", "len(r['targets']) == 8")
old = "    plan = read(WORK / 'execution-plan.json')"
assert text.count(old) == 1
text = text.replace(old, """    query = WORK / 'placement'
    assert read(query / 'placement-query-bridge.json')['accepted']
    assert read(query / 'retirement-run.json')['accepted'] and read(query / 'native-retirement.json')['pass']
    assert read(WORK / 'actual-placement-review.json')['accepted']
    assert read(WORK / 'unintended-mandatory-stop.json')['live_members_absent_after_recovery']
    assert read(WORK / 'unintended-mandatory-stop.json')['accepted_as_qualification'] is False
    plan = read(WORK / 'execution-plan.json')""")
start = text.index('    useful_names = ')
end = text.index('    blobs = {}', start)
text = text[:start] + """    useful_names = {'results.json', 'actual-results-review.json', 'actual-placement-review.json',
                    'placement-analysis.json', 'controller-review.json', 'analyzer-source-review.json',
                    'execution-plan.json', 'native-stage-preparation.json', 'source-review.json',
                    'source-review-final.json', 'host-source-review.json', 'compiled-review.json',
                    'compiler.json', 'compiler-controller-review.json', 'summary-source-review.json',
                    'stage-source-review.json', 'native-controller-source-review.json',
                    'local-checks.json', 'final-retirement-source-review.json',
                    'unintended-mandatory-stop.json', 'prior-source-disposition.json'}
    useful = [p for p in files() if (p.suffix in ('.py', '.ps1', '.cs', '.slang') and
              not {'native-temp', 'owner-temp'}.intersection(p.relative_to(WORK).parts)) or p.name in useful_names]
""" + text[end:]
start = text.index('        message = ')
end = text.index("        commit = git('commit-tree'", start)
text = text[:start] + """        message = ('Preserve reviewed native RHS layer observations and findings (#44)\\n\\n'
                   'Exact diagnostic host/source wiring at24a8e14 retains original calculations, joins, '
                   'publication, refusal/recovery and unchanged8MiB/governors. '
                   'Driver query binds seven clock sites to completeISA/ELF. '
                   'Three original native controls pass with actual input snapshots; '
                   '32 original timestamp-buffer digests remain exact. '
                   'Explicit Linux build/10authority controls and three-platform compile/nonrender CI pass. '
                   'Interrupted default Mandatory output is invalidated. '
                   'Layer intervals are inclusive perturbed rawticks, not exclusive opcode/barrier cost or ordinary benefit. '
                   'Scientific rendering/CPU integration/release qualification remain open. '
                   'Complete raw/build evidence remains local in ignored attestations, without remote raw-backup claim.\\n')
""" + text[end:]
ast.parse(text)
target = WORK / 'preserve_batch.py'
assert not target.exists()
target.write_text(text)
(WORK / 'preservation-preparation.json').write_text(json.dumps({'basis':str(OLD.relative_to(ROOT)),
    'basis_seal':{'bytes':OLD.stat().st_size,'sha256':hashlib.sha256(OLD.read_bytes()).hexdigest()},
    'current_seal':{'bytes':target.stat().st_size,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()},
    'scope':'Source adaptation only. Four-owner recovery, reviewed actual results and publication/final retirement prerequisites retained; no archive or removal yet.'},indent=2)+'\n')
print(json.dumps({'bytes':target.stat().st_size,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()}))
