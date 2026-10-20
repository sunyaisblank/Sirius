import hashlib
import json
import re
import struct
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
OLD = ROOT / 'attestations/qualification/2c8ef086/native-wave32-driver-insight/files/out/native-shader-driver-introspection/transport-default-query'
NEW = ROOT / 'out/native-coefficient-sharing-driver/transport-default-query'


def read(folder):
    data = (folder / 'shader/shader-info-1.bin').read_bytes()
    assert data[:6] == b'\x7fELF\x02\x01'
    header = struct.unpack_from('<16sHHIQQQIHHHHHH', data)
    rows = [struct.unpack_from('<IIQQQQIIQQ', data, header[6] + i * header[11])
            for i in range(header[12])]
    names_row = rows[header[13]]
    names = data[names_row[4]:names_row[4] + names_row[5]]
    sections = {names[row[0]:].split(b'\0', 1)[0].decode(): data[row[4]:row[4] + row[5]]
                for row in rows}
    text = sections['.text']
    isa = (folder / 'shader/shader-info-2.isa.txt').read_text()
    words = re.findall(r';\s+((?:[0-9a-f]{8})(?: [0-9a-f]{8})*)\s*$', isa, re.M)
    encoded = b''.join(struct.pack('<I', int(word, 16)) for row in words for word in row.split())
    assert text.startswith(encoded)
    padding = text[len(encoded):]
    assert not padding or padding == struct.pack('<I', 0xbf9f0000) * (len(padding) // 4)
    llvm = sections['.AMDGPU.comment.llvmir'].decode()
    start = llvm.index('\n2058:')
    end = llvm.index('\nCacheCoefficient.exit:', start)
    block = llvm[start:end]
    result = json.loads((folder / 'result.json').read_text(encoding='utf-8-sig'))
    metadata = json.loads((folder / 'decoded-driver-object.json').read_text())
    stage = metadata['PAL_metadata']['amdpal.pipelines'][0]['.hardware_stages']['.cs']
    assert result['dispatches'] == 0
    assert stage['.wavefront_size'] == 64 and stage['.threadgroup_dimensions'] == [64, 1, 1]
    return {'module_sha256': result['shader_sha256'], 'binary_sha256': hashlib.sha256(data).hexdigest(),
            'text_bytes': len(text), 'ISA_encoded_bytes': len(encoded), 'trailing_padding_bytes': len(padding),
            'complete_ISA_encoding_equals_text_prefix': True, 'text_sha256': hashlib.sha256(text).hexdigest(),
            'LLVM_bytes': len(sections['.AMDGPU.comment.llvmir']), 'LLVM_lines': len(llvm.splitlines()),
            'coefficient_block_lines': len(block.splitlines()),
            'coefficient_reciprocal_static_sites': block.count(' = fdiv '),
            'total_reciprocal_static_sites': llvm.count(' = fdiv '),
            'vgprs': stage['.vgpr_count'], 'sgprs': stage['.sgpr_count'],
            'LDS_bytes': stage['.lds_size'], 'scratch_bytes': stage['.scratch_memory_size'],
            'modules': result['loaded_modules']}


old, new = read(OLD), read(NEW)
assert old['module_sha256'] == '2caabcd76c28c5340c21ce5777a58b62b31bb2ed61c1218d15bb10edf0304752'
assert new['module_sha256'] == 'e634c8a246f3305a5d097e83583456a572e8fd6b1902ba8c402c9f1d1fc9b2e5'
assert [(row['name'], row['bytes'], row['sha256']) for row in old['modules']] == [(row['name'], row['bytes'], row['sha256']) for row in new['modules']]
assert (old['coefficient_reciprocal_static_sites'], new['coefficient_reciprocal_static_sites']) == (6, 2)
assert old['LDS_bytes'] == new['LDS_bytes'] and old['scratch_bytes'] == new['scratch_bytes'] == 0
report = {'schema': 'sirius-native-coefficient-sharing-actual-driver-cost-v1', 'passed': True,
          'baseline': old, 'candidate': new,
          'text_reduction_percent': 100 * (old['text_bytes'] - new['text_bytes']) / old['text_bytes'],
          'decision': 'Distinct material static code reduction survives driver compilation; one original matched-producer/runtime gate is justified after source/PE admission.',
          'countercost': 'VGPR count increases from55 to63. No dynamic occupancy or performance inference; required original benefit gates may reject this candidate.',
          'limits': 'Saved baseline and one fresh candidate diagnostic Vulkan1.2/default64/no optional features query on exact loader/provider. No dispatches; this is not application Vulkan1.3 pipeline identity, correctness or benefit evidence.'}
(WORK / 'driver-cost.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({key: value for key, value in report.items() if key not in ('baseline', 'candidate')}))
