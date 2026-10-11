"""Inspect complete existing driver payloads, without compiler or native launch."""
import sys
sys.dont_write_bytecode = True
from collections import Counter
from pathlib import Path
import hashlib, json, re, struct

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
PRIOR = ROOT / 'attestations/native-vulkan/24a8e14/native-rhs-layer-clock-execution/files/placement'

def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

def inspect(directory):
    result = json.loads((directory / 'clock-placement.json').read_text())
    observation, = result['observations']
    for query in observation['queries']:
        assert query['captured'] and query['available']
        assert query['size_query_result'] == query['complete_query_result'] == 0
        assert seal(directory / 'query-payloads' / query['file']) == {
            'bytes': query['returned_bytes'], 'sha256': query['sha256']}
    data = (directory / 'query-payloads/transport.info1').read_bytes()
    assert data[:6] == b'\x7fELF\x02\x01'
    offset = struct.unpack_from('<Q', data, 40)[0]
    stride, count, names_index = struct.unpack_from('<HHH', data, 58)
    assert stride == 64 and names_index < count
    headers = [struct.unpack_from('<IIQQQQIIQQ', data, offset + i * stride)
               for i in range(count)]
    header = headers[names_index]
    names = data[header[4]:header[4] + header[5]]
    sections = {}
    for header in headers:
        name = names[header[0]:names.index(0, header[0])].decode('utf-8')
        assert name not in sections
        sections[name] = data[header[4]:header[4] + header[5]]
    raw = (directory / 'query-payloads/transport.info2').read_bytes()
    assert raw[-1] == 0
    instructions, labels, printed = [], {}, bytearray()
    for line_number, line in enumerate(raw[:-1].decode('utf-8').splitlines(), 1):
        label = re.fullmatch(r'(_L[0-9a-f]+):', line.strip())
        if label:
            labels[label[1]] = len(printed)
        match = re.search(r';\s*((?:[0-9a-f]{8}\s*)+)$', line)
        if not match:
            continue
        words = [int(word, 16) for word in match[1].split()]
        instructions.append({'pc': len(printed), 'line': line_number,
                             'text': line.split(';')[0].strip(), 'words': words})
        printed.extend(b''.join(struct.pack('<I', word) for word in words))
    assert bytes(printed) == sections['.text'][:len(printed)]
    branches = []
    for inst in instructions:
        match = re.fullmatch(r'(s_(?:cbranch_\w+|branch))\s+(_L[0-9a-f]+)', inst['text'])
        if not match:
            continue
        assert len(inst['words']) == 1
        displacement = inst['words'][0] & 0xffff
        if displacement & 0x8000:
            displacement -= 0x10000
        target = inst['pc'] + 4 + 4 * displacement
        assert target == labels[match[2]]
        branches.append({**inst, 'target': target, 'backward': target <= inst['pc']})
    reciprocals = []
    for inst in instructions:
        if inst['text'].startswith('v_rcp_f32'):
            enclosing = sorted((branch for branch in branches if branch['backward'] and
                                branch['target'] <= inst['pc'] <= branch['pc']),
                               key=lambda branch: branch['pc'] - branch['target'])
            reciprocals.append({**inst, 'enclosing_backward_edges': enclosing})
    llvm = sections['.AMDGPU.comment.llvmir'].decode('utf-8')
    return {'module': {k: observation[k] for k in ('module_bytes', 'module_sha256')},
            'statistics': observation['queries'][0]['statistics'],
            'text_bytes': len(sections['.text']), 'printed_text_bytes': len(printed),
            'all_printed_bytes_match_elf_prefix': True,
            'instruction_count': len(instructions),
            'opcode_counts': dict(sorted(Counter(i['text'].split()[0] for i in instructions).items())),
            'validated_direct_branches': len(branches), 'reciprocal_sites': reciprocals,
            'llvm_bytes': len(sections['.AMDGPU.comment.llvmir']),
            'wave64_attribute': '+wavefrontsize64' in llvm,
            'branch_edges': branches, 'instructions': instructions}

assert WORK == ROOT / 'out/native-bound-selection'
bridge = json.loads((WORK / 'placement/placement-query-bridge.json').read_text())
assert bridge['accepted'] and bridge['returncode'] == 0
baseline, candidate = inspect(PRIOR), inspect(WORK / 'placement')
assert baseline['module']['module_sha256'] == 'e3a1589d852a4254312bb0e9e4b0343ab6aea32855c02c4f0f93b9fe1e18e5a3'
assert candidate['module']['module_sha256'] == '0cb8acec4646038f184d2553d8def1f3f5cb538ab0c2fcf73f7ab82c001a5d37'
record = {'accepted': True, 'source_revision': bridge['source_revision'],
          'baseline': baseline, 'candidate': candidate,
          'scope': 'Whole emitted static code and resource comparison for perturbed '
          'Transport observer pipeline. Direct branch encodings checked against ELF. '
          'Backward-edge lexical ranges are not complete CFGs or exclusive costs. '
          'Zero dispatch; no numerical identity or complete-path benefit verdict.'}
path = WORK / 'driver-inspection.json'
assert not path.exists()
path.write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({name: {key: record[name][key] for key in (
    'statistics', 'text_bytes', 'printed_text_bytes', 'instruction_count',
    'validated_direct_branches', 'wave64_attribute')} for name in ('baseline', 'candidate')}))
