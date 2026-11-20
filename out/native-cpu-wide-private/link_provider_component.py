"""Link the provider additions against exactly one copy of the AX math cores."""

from pathlib import Path
import os
import subprocess
import sys

from prepare_provider import AX, R, fingerprint_tree
from replay_native import EMIT, SEALED
from resume_raw import check_space, pin, read_json, write_json
from run_value_abi import expect_pin


def symbols(path, *flags):
    return {line.split()[-1] for line in subprocess.check_output(
        ['/usr/bin/nm', *flags, str(path)], text=True).splitlines()}


def main():
    scratch, = (Path(arg).resolve() for arg in sys.argv[1:])
    if scratch == R / 'out' or not scratch.is_relative_to(R / 'out'):
        raise RuntimeError('Use the working directory for this provider build under out/')
    check_space()
    prepared = read_json(scratch / 'provider-inputs.json')
    if fingerprint_tree(scratch / 'mesa') != prepared['candidate']:
        raise RuntimeError('The prepared Mesa source identity changed')
    adapter = scratch / 'admitted_rows.o'
    expect_pin(adapter, prepared['adapter'])
    proof = read_json(SEALED)['results']
    replay = read_json(scratch / 'native-replay.json')['results']
    ax = read_json(AX)
    ax_objects = [R / row['path'] for row in ax['native_objects']]
    for path, record in zip(ax_objects, ax['native_objects'], strict=True):
        expect_pin(path, record)

    # The sealed SIMD16 diagnostic component also links these same nine cores.
    # At the provider boundary AX already owns them, so link only the wide
    # adapter and stages, with an explicit check of their remaining references.
    core = R / ax['closed_dp_math_object']['path']
    expect_pin(core, ax['closed_dp_math_object'])
    names = {name for name in proof['link']['symbol_audit']['defined_symbols']
             if name.startswith('sirius_local_value_')}
    if len(names) != 9 or not names <= symbols(ax_objects[-1], '-g', '--defined-only'):
        raise RuntimeError('The original provider component lacks an exact math core')
    if {name for name in symbols(core, '-g', '--defined-only') if name.startswith('sirius_local_value_')} != names:
        raise RuntimeError('The shared math object defines a different core set')

    parts = [adapter]
    for action in (*EMIT, 'adapter-object'):
        record = replay[action]['output']
        path = R / record['path']
        expect_pin(path, proof[action]['output'])
        parts.append(path)
    output = R / 'bin/linux-gcc/native-cpu-wide-private/wide-provider.o'
    pending = output.with_suffix('.pending.o')
    if pending.exists():
        raise RuntimeError('An unfinished component link exists')
    if output.exists():
        expect_pin(output, prepared['adapter'])
    argv = ['/usr/bin/ld', '-r', *map(str, parts), '-o', str(pending)]
    subprocess.run(argv, check=True)
    undefined = sorted(symbols(pending, '-u'))
    expected = names | set(proof['link']['symbol_audit']['external_symbols'])
    if set(undefined) != expected:
        raise RuntimeError(f'Unexpected provider link dependencies: {sorted(set(undefined) ^ expected)}')
    os.replace(pending, output)
    write_json(scratch / 'provider-component.json', {
        'source': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'output': {'path': str(output.relative_to(R)), **pin(output)},
        'parts': [{'path': str(path.relative_to(R)), **pin(path)} for path in parts],
        'shared_math': ax['closed_dp_math_object'], 'shared_math_owner': ax['native_objects'][-1],
        'undefined': undefined, 'argv': argv,
    })


if __name__ == '__main__':
    main()
