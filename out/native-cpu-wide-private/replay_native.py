"""Recreate the sealed SIMD16 objects for a provider build in this checkout."""

from pathlib import Path
import subprocess
import sys

from recover_sdk_ranges import recover
from resume_raw import (BASELINE, R, acquire_ispc, check_space, emit, generated_sources,
                        read_json, recover_sources, write_json)
from run_value_abi import expect_pin


HERE = Path(__file__).resolve().parent
SEALED = R / 'attestations/qualification/19fba3ef/native-cpu-simd16-native/MANIFEST.json'
EMIT = ('emit-transport', 'emit-endpoint', 'emit-dopri_phase')


def main():
    scratch, = (Path(arg).resolve() for arg in sys.argv[1:])
    if scratch == R / 'out' or not scratch.is_relative_to(R / 'out') or not scratch.is_dir():
        raise RuntimeError('Pass an existing per-task directory under out/')
    if Path.cwd() != R or subprocess.check_output(['git', 'status', '--porcelain']).strip():
        raise RuntimeError('Run from the clean Sirius checkout')
    check_space()
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    baseline = read_json(BASELINE)
    inputs = baseline['replay']['source_inputs']
    recover_sources(inputs)
    acquire_ispc(scratch, inputs)
    contract = generated_sources(inputs, baseline['restored_generated_sources'])
    raw = {}
    emit(scratch, revision, baseline, contract, raw, EMIT)
    write_json(scratch / 'raw-replay.json', raw)
    if not (scratch / 'sdk').is_dir():
        recover(scratch)
    for script, args in (
        ('run_value_abi.py', ('value-replay.json',)),
        ('run_context.py', ('context-replay.json',)),
        ('run_native_context.py', ('context-replay.json', 'native-replay.json')),
    ):
        subprocess.run([sys.executable, '-B', str(HERE / script), str(scratch),
                        *(str(scratch / arg) for arg in args)], check=True)
    sealed = read_json(SEALED)['results']
    replay = read_json(scratch / 'native-replay.json')['results']
    for name in (*EMIT, 'adapter-object', 'link'):
        path = R / replay[name]['output']['path']
        expect_pin(path, sealed[name]['output'])


if __name__ == '__main__':
    main()
