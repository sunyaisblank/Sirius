"""Inspect the existing compiler outputs without rerunning any compiler."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, re, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
PRIOR = ROOT / 'attestations/native-vulkan/24a8e14/native-rhs-layer-clock-execution/files'

def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

def parse(text):
    functions = re.findall(r'(%\S+) = OpFunction\b.*?OpFunctionEnd', text, re.S)
    bodies = re.findall(r'%\S+ = OpFunction\b.*?OpFunctionEnd', text, re.S)
    assert len(functions) == len(bodies) == len(set(functions))
    definitions = dict(re.findall(r'^\s*(%\S+) = (.*)$', text, re.M))
    names = dict(re.findall(r'OpName (%\S+) "([^"]+)"', text))
    function_names = {symbol: names[symbol] for symbol in functions}
    cache = {}
    def global_value(symbol):
        if symbol in function_names:
            return 'FUNCTION:' + function_names[symbol]
        if symbol not in cache:
            value = definitions[symbol]
            cache[symbol] = re.sub(r'%\S+', lambda m: global_value(m[0]), value)
        return cache[symbol]
    canonical = {}
    for symbol, body in zip(functions, bodies):
        name = names[symbol]
        local = {s: 'LOCAL' + str(i) for i, s in enumerate(
            re.findall(r'^\s*(%\S+) = ', body, re.M))}
        local.pop(symbol)
        normalized = re.sub(r'%\S+', lambda m: local[m[0]] if m[0] in local
                            else global_value(m[0]), body)
        normalized = '\n'.join(line.strip() for line in normalized.splitlines())
        canonical.setdefault(name, []).append(normalized)
    return canonical

assert WORK == ROOT / 'out/native-bound-selection'
receipt = json.loads((WORK / 'compiler.json').read_text())
assert not receipt['accepted'] and receipt['errors'] == ['AssertionError()']
assert len(receipt['runs']) == 4 and all(r['accepted'] for r in receipt['runs'])
assert not receipt['signals']
for path, expected in receipt['input_seals'].items():
    assert seal(Path(path)) == expected
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == receipt['baseline_revision']
assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
baseline_path = PRIOR / 'phase/retained_transport.spvasm'
current_path = WORK / 'phase/retained_transport.spvasm'
before, after = baseline_path.read_text(), current_path.read_text()
assert after.count('OpReadClockKHR') == before.count('OpReadClockKHR') == 2
assert after.count('OpControlBarrier') == before.count('OpControlBarrier') == 9
assert after.count('NoContraction') == before.count('NoContraction')
assert 'OpCapability Int64' not in after and 'OpTypeInt 64' not in after
assert 'OpTypeFloat 64' not in after
assert not re.search(r'FPFastMathMode|RelaxedPrecision|OpFmaKHR', after)
decorated = set(re.findall(r'OpDecorate (%\S+) NoContraction', after))
helpers = {}
for name, expected in (('RPUpAdd', 'OpFAdd'), ('RPUpMultiply', 'OpFMul')):
    body = re.search(r'%' + name + r' = OpFunction\b.*?OpFunctionEnd', after, re.S)[0]
    assert not re.search(r'OpBranchConditional|OpSelectionMerge|OpLoopMerge', body)
    operations = re.findall(r'(%\S+) = (OpF(?:Add|Sub|Mul|Div)) ', body)
    assert len(operations) == 1 and operations[0][1] == expected
    assert operations[0][0] in decorated
    helpers[name] = {'branches': 0, 'ordered_operation': expected,
                     'NoContraction': True, 'body': body}
old, new = parse(before), parse(after)
changed = [name for name in sorted(old.keys() & new.keys()) if old[name] != new[name]]
removed, added = sorted(old.keys() - new.keys()), sorted(new.keys() - old.keys())
assert changed == ['RPUpAdd', 'RPUpMultiply']
assert not removed and not added
record = {'accepted': True, 'scope': 'Pure inspection of actual four successful tool outputs; '
          'original postcheck assumed calls which Slang inlined. Raw failed postcheck preserved. '
          'No repeat compile, device/dispatch, numerical or benefit verdict.',
          'source_revision': receipt['baseline_revision'],
          'compiler_receipt': seal(WORK / 'compiler.json'),
          'baseline_assembly': seal(baseline_path), 'assembly': seal(current_path),
          'spv': seal(WORK / 'phase/retained_transport.spv'),
          'helpers': helpers, 'changed_bodies': changed,
          'removed_unused_functions': removed, 'unchanged_function_bodies':
          sum(len(old[n]) for n in old.keys() & new.keys() if n not in changed),
          'barriers': 9, 'clock_instructions': 2,
          'NoContraction_count': after.count('NoContraction'),
          'input_seals_unchanged': len(receipt['input_seals'])}
target = WORK / 'compile-admission.json'
assert not target.exists()
target.write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({k: v for k, v in record.items() if k != 'helpers'}))
