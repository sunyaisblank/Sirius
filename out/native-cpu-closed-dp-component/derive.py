"""Derive only AU's admitted DP v4 callback; never compile or patch Mesa."""
from pathlib import Path
import hashlib
import json

R = Path.cwd()
W = Path(__file__).resolve().parent
AU = R / 'out/native-cpu-shared-value-component'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def main():
    proof = json.loads((AU / 'derivation.json').read_text())
    original = AU / 'rows_adapter.cpp'
    assert pin(original) == {k: proof['derived'][k] for k in ('bytes', 'sha256')}
    source = original.read_text()
    before, after = 'sirius_scope_dp_packet(', 'sirius_closed_dp_packet('
    assert source.count(before) == 2 and after not in source
    candidate = source.replace(before, after)
    assert candidate.replace(after, before) == source
    assert candidate[:source.index('extern "C" __attribute__((visibility("hidden"))) void sirius_scope_transport_packet')]==source[:source.index('extern "C" __attribute__((visibility("hidden"))) void sirius_scope_transport_packet')]
    target = W / 'rows_adapter.cpp'
    with target.open('x') as stream:
        stream.write(candidate)
    command = list(proof['compile_argv'])
    command[command.index('-c') + 1] = str(target)
    command[command.index('-o') + 1] = str(R / 'bin/linux-gcc/native-cpu-closed-dp-component/adapter-only.o')
    record = dict(prepared_only=True, compiler_run=False,
        source_inputs={str(path.relative_to(R)): pin(path) for path in (original, AU / 'derivation.json')},
        derived=dict(path=str(target.relative_to(R)), **pin(target)), compile_argv=command,
        source_inverse_exact=True, public_wrapper_prefix_byte_preserved=True,
        private_change='Only the admitted DP v4 declaration/call changes from sirius_scope_dp_packet to sirius_closed_dp_packet. AU Transport/Endpoint and all public/scalar/AVX2/old-v4 branches and guards remain exact.',
        scope='Source only. The closed DP body requires the genuine positive current job, original pure preflight and complete saved/masked/RNE/DAZ-FTZ-off/restored MXCSR scope. No numerical acceptance transfer; original C8/all32 is the next real gate. No ABBA/frame authorization until remaining-stage closure is justified.')
    with (W / 'derivation.json').open('x') as stream:
        json.dump(record, stream, indent=2)
        stream.write('\n')


if __name__ == '__main__':
    main()
