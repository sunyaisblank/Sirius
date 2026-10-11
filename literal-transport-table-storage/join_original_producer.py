"""Relink exact current Ninja inputs and join the accepted b81 product bytes."""
from pathlib import Path
import hashlib, json, shlex, subprocess, sys
sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[1]
BUILD = ROOT / 'bin/linux-gcc'
ARCHIVE = ROOT / 'attestations/review/b81e061/rejected-upmultiply-restoration-source-checks'

def seal(path):
    data = Path(path).read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

def main():
    receipt_path = ARCHIVE / 'diagnostic/restoration/linux-readonly-payloads.json'
    receipt = json.loads(receipt_path.read_text())
    assert receipt['source_revision'] == 'b81e061fa081cdaed37eca3add950bd254bf296d' and receipt['pass_']
    original = BUILD / 'tests/backend/sirius_backend_tests'
    expected = receipt['actual_consumers']['bin/linux-gcc/tests/backend/sirius_backend_tests']
    expected = {key: expected[key] for key in ('bytes', 'sha256')}
    assert seal(original) == expected
    commands = subprocess.check_output(['ninja', '-C', str(BUILD), '-t', 'commands', 'sirius_backend_tests'], text=True)
    lines = [line for line in commands.splitlines() if ' -o tests/backend/sirius_backend_tests ' in line]
    assert len(lines) == 1
    words = shlex.split(lines[0]); first = words.index('&&')+1; end = words.index('&&', first)
    argv = words[first:end]
    assert argv[0] == '/usr/bin/g++-14' and argv.count('-o') == 1
    at = argv.index('-o')+1; assert argv[at] == 'tests/backend/sirius_backend_tests'
    inputs = [Path(arg) if Path(arg).is_absolute() else BUILD / arg for i,arg in enumerate(argv) if i != at and (arg.endswith(('.a','.o','.so')) or i == 0)]
    before = {str(p.resolve(strict=True)): seal(p) for p in inputs}
    destination = WORK / 'canonical-reference-backend-tests'; assert not destination.exists()
    argv[at] = str(destination)
    with (WORK/'reference-link.stdout').open('wb') as out, (WORK/'reference-link.stderr').open('wb') as err:
        result = subprocess.run(argv, cwd=BUILD, stdout=out, stderr=err)
    after = {name: seal(name) for name in before}
    report = {'scope':'Exact relink to accepted source/build product using current original Ninja object/archive/link command; no independent recompile of reused objects or scientific execution', 'argv':argv, 'cwd':str(BUILD), 'exit':result.returncode, 'inputs_before':before, 'inputs_after':after, 'receipt':{'path':str(receipt_path), **seal(receipt_path)}, 'accepted_product':{'path':str(original), **expected}, 'relinked_product':{'path':str(destination), **seal(destination)} if destination.exists() else None, 'pass':result.returncode==0 and before==after and destination.exists() and seal(destination)==expected and seal(original)==expected}
    (WORK/'original-producer-join.json').write_text(json.dumps(report,indent=2)+'\n')
    assert report['pass'], report
    print(json.dumps({'pass':True,'original_relink_byte_exact_accepted_b81':True,'link_inputs':len(before)}))

if __name__ == '__main__':
    main()
