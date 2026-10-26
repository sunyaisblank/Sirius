"""Check retained scalar shader code and the two unaffected pure predicates."""
from pathlib import Path
import json
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
sys.path.insert(0, str(R / 'out/native-cpu-field-compatibility'))
from elf_proof import Elf, pin


def functions(elf):
    return {s['name']: s for s in elf.symbols
            if s['section'] and s['size'] and s['info'] & 15 == 2}


def relocations(elf, symbol):
    return [(offset - symbol['value'], kind, addend, target['name'])
            for section, offset, kind, addend, target in elf.relocations
            if section == symbol['section'] and 0 <= offset - symbol['value'] < symbol['size']]


rows = []
for old, new, excluded in (
    ('bin/linux-gcc/native-cpu-scalar-zero/objects/endpoint.o',
     'bin/linux-gcc/native-cpu-endpoint-header-domain/endpoint.o',
     {'sirius_native_endpoint_try_packet'}),
    ('bin/linux-gcc/native-cpu-dispatch-admission/admission.o',
     'bin/linux-gcc/native-cpu-endpoint-header-domain/admission.o',
     {'sirius_native_endpoint_admit_packet'})):
    a, b = Elf(R / old), Elf(R / new)
    af, bf = functions(a), functions(b)
    assert set(af) == set(bf)
    same = sorted(set(af) - excluded)
    for name in same:
        assert a.function_bytes(af[name]) == b.function_bytes(bf[name]), name
        assert relocations(a, af[name]) == relocations(b, bf[name]), name
    assembly = subprocess.check_output(['objdump', '-d', '-M', 'intel', '--no-show-raw-insn',
                                        str(R / new)], text=True)
    mnemonics = [line.split(':', 1)[1].strip().split()[0] for line in assembly.splitlines()
                 if ':' in line and line.split(':', 1)[0].strip().isalnum()
                 and line.split(':', 1)[1].strip()]
    assert not any(word.startswith('v') for word in mnemonics), new
    rows.append(dict(original=dict(path=old, **pin(R / old)),
                     candidate=dict(path=new, **pin(R / new)),
                     unchanged_function_count=len(same), unchanged_functions=same,
                     function_bytes_and_relocation_targets_identical=True,
                     changed_functions=sorted(excluded), optional_ISA_instructions=0))
with (W / 'code-proof.json').open('x') as stream:
    json.dump(dict(passed=True, objects=rows,
        scope='166 original scalar shader/helper functions and Transport/DopriPhase pure predicates retain exact bytes and relocation tuples. Only two Endpoint routing wrappers change. Baseline ISA only. No numerical, provider, timing or release acceptance.'), stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(passed=True, unchanged_functions=sum(row['unchanged_function_count'] for row in rows))))
