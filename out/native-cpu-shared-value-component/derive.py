"""Prepare only the private admitted callbacks for the coupled value candidate."""
from pathlib import Path
import hashlib
import json

R = Path.cwd()
W = Path(__file__).resolve().parent
old = R / 'out/native-cpu-dispatch-admission-token'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,
                    sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


proof = json.loads((old / 'derivation.json').read_text())
source = (old / 'rows_adapter.cpp').read_text()
expected = proof['derived']['rows_adapter.cpp'] if 'rows_adapter.cpp' in proof['derived'] else proof['derived']['out/native-cpu-dispatch-admission-token/rows_adapter.cpp']
assert pin(old / 'rows_adapter.cpp') == expected
boundary = source.index('#include "admitted_bridge.h"')
public = source[:boundary]
private = source[boundary:]
declarations = []
for stage, name in (('transport', 'transport'), ('endpoint', 'endpoint'), ('dopri_phase', 'dp')):
    before = 'sirius_simd_v4_' + stage + '_packet('
    after = 'sirius_scope_' + name + '_packet('
    assert private.count(before) == 1
    private = private.replace(before, after, 1)
    declarations.append('extern "C" __attribute__((visibility("hidden"))) void ' +
        after + 'uint32_t *input, uint32_t input_words, uint32_t *output, uint32_t output_words, uint32_t first_x, uint32_t groups);')
private = private.replace('#include "admitted_bridge.h"', '#include "../native-cpu-dispatch-admission-token/admitted_bridge.h"', 1)
private = private.replace('#include "native_packet.hpp"', '#include "../native-cpu-dispatch-admission-token/native_packet.hpp"', 1)
candidate = public + '\n'.join(declarations) + '\n' + private
assert candidate[:boundary] == public
target = W / 'rows_adapter.cpp'
with target.open('x') as stream:
    stream.write(candidate)
command, = [r['argv'] for r in proof['prospective_compile_commands'] if r['role'] == 'rows_adapter']
command = list(command)
command[command.index('-c') + 1] = str(target)
command[command.index('-o') + 1] = str(R / 'bin/linux-gcc/native-cpu-shared-value-component/rows_adapter.o')
record = dict(prepared_only=True, compiler_run=False,
    source_inputs={str((old / n).relative_to(R)): pin(old / n)
                   for n in ('derivation.json', 'rows_adapter.cpp', 'admitted_bridge.h', 'native_packet.hpp')},
    derived=dict(path=str(target.relative_to(R)), **pin(target)), compile_argv=command,
    public_wrapper_prefix_byte_preserved=True,
    private_change='Only three admitted v4 raw callbacks change to their separately named strict value modules. All existing admission/ISA/FP/grid/range/alias/row guards and scalar/AVX2 fallbacks remain exact.',
    scope='Source preparation only; no runtime scope or numerical compatibility/adoption inferred. Required genuine job flag/full-MXCSR scope and coupled native3+shared9core link must be proved before qualification.')
with (W / 'derivation.json').open('x') as stream:
    json.dump(record, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(prepared=True, public_prefix_exact=True, changed_private_callbacks=3)))
