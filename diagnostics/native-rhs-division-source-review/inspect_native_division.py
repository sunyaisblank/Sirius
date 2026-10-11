"""Inspect sealed existing IR/ISA only; no compilation or shader execution."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib
import json
import re
import struct
import subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
ARCHIVE = ROOT / 'attestations/native-vulkan/24a8e14/native-rhs-layer-clock-execution'


def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def elf_text(data):
    assert data[:6] == b'\x7fELF\x02\x01'
    offset = struct.unpack_from('<Q', data, 40)[0]
    stride, count, names = struct.unpack_from('<HHH', data, 58)
    assert stride == 64 and names < count
    headers = [struct.unpack_from('<IIQQQQIIQQ', data, offset + i * stride)
               for i in range(count)]
    h = headers[names]
    strings = data[h[4]:h[4] + h[5]]
    found = []
    for h in headers:
        end = strings.index(0, h[0])
        if strings[h[0]:end] == b'.text':
            found.append(data[h[4]:h[4] + h[5]])
    assert len(found) == 1
    return found[0]


assert WORK == ROOT / 'out/native-rhs-division-source-review'
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == \
    '2053ef79c90752d37b0075b8137cfb9b571b6820'
assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
manifest = json.loads((ARCHIVE / 'manifest.json').read_text())
assert seal(ARCHIVE / 'manifest.json')['sha256'] == \
    '71e2994b331a5a1e8aef3e3e5e472e19b1ce72e9b9c2aaed7430efff08510d2d'
paths = ('phase/retained_transport.spvasm',
         'placement/query-payloads/transport.info1',
         'placement/query-payloads/transport.info2')
inputs = {}
for name in paths:
    path = ARCHIVE / 'files' / name
    row, = [r for r in manifest['records'] if r['recovery'] == str(path.relative_to(ROOT))]
    inputs[name] = seal(path)
    assert inputs[name] == {k: row[k] for k in ('bytes', 'sha256')}
ir = (ARCHIVE / 'files' / paths[0]).read_text()
functions = {}
for name in ('EvaluateOne', 'RTDivide', 'RTDivide_0', 'RPDivideSeed', 'RPRatioUpper',
             'RPFloorExponent', 'RPMantissa'):
    match = re.search(r'^\s*%' + re.escape(name) + r' = OpFunction\b.*?OpFunctionEnd',
                      ir, re.M | re.S)
    assert match
    text = match[0]
    functions[name] = {'lines': len(text.splitlines()),
                       'calls': [line.strip() for line in text.splitlines()
                                 if 'OpFunctionCall' in line],
                       'source': text}
assert any('%RTDivide_0 ' in s for s in functions['EvaluateOne']['calls'])
code = elf_text((ARCHIVE / 'files' / paths[1]).read_bytes())
raw = (ARCHIVE / 'files' / paths[2]).read_bytes()
assert raw[-1] == 0
printed = bytearray()
instructions = []
labels = {}
for line_number, line in enumerate(raw[:-1].decode('utf-8').splitlines(), 1):
    label = re.fullmatch(r'(_L[0-9a-f]+):', line.strip())
    if label:
        labels[label[1]] = len(printed)
    match = re.search(r';\s*((?:[0-9a-f]{8}\s*)+)$', line)
    if not match:
        continue
    words = [int(word, 16) for word in match[1].split()]
    instructions.append({'line': line_number, 'pc': len(printed),
                         'text': line.split(';')[0].strip(), 'words': words})
    printed.extend(b''.join(struct.pack('<I', word) for word in words))
assert bytes(printed) == code[:len(printed)]
assert len(printed) == 237116 and len(code) == 237568
branches = []
for inst in instructions:
    match = re.fullmatch(r'(s_(?:cbranch_\w+|branch))\s+(_L[0-9a-f]+)', inst['text'])
    if not match:
        continue
    assert len(inst['words']) == 1 and match[2] in labels
    displacement = inst['words'][0] & 0xffff
    if displacement & 0x8000:
        displacement -= 0x10000
    target = inst['pc'] + 4 + 4 * displacement
    assert target == labels[match[2]]
    branches.append({**inst, 'target': target, 'backward': target <= inst['pc']})
reciprocals = []
for index, inst in enumerate(instructions):
    if not inst['text'].startswith('v_rcp_f32'):
        continue
    enclosing = [b for b in branches if b['backward'] and
                 b['target'] <= inst['pc'] <= b['pc']]
    enclosing.sort(key=lambda b: b['pc'] - b['target'])
    reciprocals.append({**inst, 'between_layer_body_cuts': 0x8220 < inst['pc'] < 0x1a02c,
                        'enclosing_direct_backward_edges': enclosing,
                        'context': instructions[max(0, index - 65):index + 28]})
assert len(reciprocals) == 9
result = {'production_source': '2053ef79c90752d37b0075b8137cfb9b571b6820',
          'diagnostic_source': manifest['source_revision'], 'inputs': inputs,
          'all_printed_bytes_match_elf_prefix': True,
          'validated_direct_branches': len(branches), 'functions': functions,
          'reciprocal_sites': reciprocals,
          'scope': 'Exact existing IR/ISA byte inspection. Enclosing direct backward edges are '
                   'lexical intervals, not a complete CFG or proof of execution/dominance. '
                   'Source attribution requires operand and branch validation. No native run, '
                   'exclusive cost, optimization benefit or qualification result.'}
target = WORK / 'native-division-inspection.json'
assert not target.exists()
target.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'direct_branches': len(branches),
                  'reciprocals': [{'pc': hex(r['pc']), 'line': r['line'],
                                   'inside_layer_body': r['between_layer_body_cuts'],
                                   'smallest_backward_edge': None if not r['enclosing_direct_backward_edges']
                                   else {k: hex(r['enclosing_direct_backward_edges'][0][k])
                                         for k in ('target', 'pc')}} for r in reciprocals]}))
