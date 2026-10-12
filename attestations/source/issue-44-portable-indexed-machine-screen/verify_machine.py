"""Join the complete offline object to the actual candidate JIT artifact."""
from pathlib import Path
import hashlib
import json
import resource
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024 * 1024, 64 * 1024 * 1024))
from elf_layout import elf

root = Path.cwd()
work = Path(__file__).resolve().parent
reference = work / 'capture-candidate/jit-10.elf'
artifact = root / 'bin/linux-gcc/portable-indexed-limb-accumulation/endpoint.o'
sha = lambda data: hashlib.sha256(data).hexdigest()
before, current = elf(reference), elf(artifact)
assert set(before[2]) == {'cs_variant', 'cs_co_variant', 'cs_co_variant.resume', 'cs_co_variant.destroy'}
assert before[2] == current[2], 'complete function layout'
assert set(before[1]) == set(current[1]), 'allocated section set'
assert before[3] == current[3], 'complete symbolic relocation records'
joins = []
for name, original in before[1].items():
    after = current[1][name]
    assert tuple(original[k] for k in ('flags', 'align', 'bytes')) == tuple(
        after[k] for k in ('flags', 'align', 'bytes')), (name, 'layout')
    assert original['data'] == after['data'], (name, 'complete bytes')
    joins.append({'section': name, 'bytes': original['bytes'], 'sha256': sha(original['data']),
                  'symbolic_relocations': len(before[3].get(name, []))})
result = {'passed': True, 'reference_sha256': sha(before[0]), 'artifact_sha256': sha(current[0]),
          'complete_functions': before[2], 'complete_allocated_sections': joins,
          'symbolic_relocations': sum(len(v) for v in before[3].values()),
          'producer_options_uniquely_recovered': False,
          'scope': 'Complete candidate JIT/offline object code, constants, function and allocated-section layout and symbolic relocation correspondence. No source-wide mapping, live PCs, arithmetic execution, exclusive cost or runtime benefit claim.'}
(work / 'machine-verification.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
