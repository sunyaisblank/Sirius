"""Resume the pinned SIMD16 ordinary IR build after removal of the preset cache."""

from collections import Counter
from pathlib import Path
import json
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile

sys.dont_write_bytecode = True
from raw_source import ISPC, PROFILES, R, pin, planned


ARCHIVES = (
    '536ac050/native-cpu-bridge-source-and-boundaries/actual-native-products.zip',
    '4332fb95/native-cpu-closed-dp-original-c8-attempt/source-and-evidence.zip',
    'd970f13e/native-cpu-closed-remaining-first-emit/source-and-ordinary-ir.zip',
)
STAGES = ('emit-endpoint', 'emit-dopri_phase')
GUARD = R / 'out/native-cpu-row-ownership/cold-frame'
BASELINE = R / 'attestations/qualification/a3b8c9f0/native-cpu-simd16-raw/MANIFEST.json'
MIN_FREE_BYTES = 100_000_000_000


def read_json(path):
    return json.loads(path.read_text())


def write_json(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def check_space():
    free = shutil.disk_usage('/').free
    if free < MIN_FREE_BYTES:
        raise RuntimeError(f'Root filesystem has only {free:,} bytes free')
    return free


def restore(source, path, expected):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as stream:
        shutil.copyfileobj(source, stream, length=1024 * 1024)
    if pin(path) != expected:
        raise RuntimeError(f'Pinned identity mismatch for {path.relative_to(R)}')


def recover_sources(inputs):
    pending = {name for name in inputs if name.endswith('.program.hpp')}
    pending.update(profile['source'] for profile in PROFILES.values())
    for name in ARCHIVES:
        with zipfile.ZipFile(R / 'attestations/qualification' / name) as archive:
            for member in sorted(pending.intersection(archive.namelist())):
                with archive.open(member) as stream:
                    restore(stream, R / member, inputs[member])
                pending.remove(member)
    if pending:
        raise RuntimeError(f'No archived source for {sorted(pending)}')


def acquire_ispc(scratch, inputs):
    receipt = read_json(R / 'out/native-cpu-simd-rows/compiler-acquisition.json')
    asset = receipt['asset']
    archive_path = scratch / asset['name']
    subprocess.run(['/usr/bin/curl', '--fail', '--location', '--silent', '--show-error',
                    '--max-time', '180', '--proto', '=https', '--output', str(archive_path),
                    asset['browser_download_url']], check=True)
    if pin(archive_path) != {'bytes': asset['size'], 'sha256': asset['digest'].removeprefix('sha256:')}:
        raise RuntimeError('ISPC release archive identity mismatch')
    members = {str(path.relative_to(R)) for path in (ISPC, ISPC.parent.parent / 'ispc.html')}
    with tarfile.open(archive_path, 'r:gz') as archive:
        for name in sorted(members):
            record = receipt['extracted_members'][name]
            member = archive.getmember(record['archive_member'])
            if not member.isfile() or member.size != inputs[name]['bytes']:
                raise RuntimeError(f'Invalid compiler archive member {member.name}')
            with archive.extractfile(member) as stream:
                restore(stream, R / name, inputs[name])
            (R / name).chmod(record['mode'])
    archive_path.unlink()


def generated_sources(inputs, expected_outputs):
    outputs, contract = planned()
    if contract['source_inputs'] != inputs or contract['outputs'] != expected_outputs:
        raise RuntimeError('Regenerated SIMD16 sources differ from the original derivation')
    for path, data in outputs.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            stream.write(data)
    return contract


def ir_inventory(path):
    ir = path.read_text()
    if re.search(r'\b(?:fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b', ir):
        raise RuntimeError(f'Fast-math operation in {path.name}')
    if re.search(r'@llvm\.(?:fma|fmuladd)\.', ir):
        raise RuntimeError(f'Fused arithmetic in {path.name}')
    calls = Counter(re.findall(r'\bcall\s+[^\n]*?@([\w.$-]+)\(', ir))
    math = {name: count for name, count in calls.items() if name.startswith('sirius_wide_pointer_')}
    context, = re.findall(r'^%v16_varying_KernelContext_0 = type .+$', ir, re.M)
    if not math or any(not name.endswith('16') for name in math):
        raise RuntimeError(f'Unexpected pointer math lane width in {path.name}')
    return {'path': str(path.relative_to(R)), **pin(path),
            'definitions': len(re.findall(r'^define ', ir, re.M)), 'context': context,
            'calls': dict(calls),
            'fp_operations': dict(Counter(re.findall(
                r'^\s+%[^=\n]+ = (fadd|fsub|fmul|fdiv|frem|fneg|fcmp|sitofp|uitofp|fptosi|fptoui|fpext|fptrunc)\b',
                ir, re.M))),
            'comparison_predicates': dict(Counter(re.findall(r'= fcmp (\w+) ', ir)))}


def emit(scratch, revision, baseline, contract, results):
    # The runner resolves its own location. Its wrapper stays in the task directory.
    wrapper = scratch / 'run_owned.py'
    shutil.copyfile(GUARD / wrapper.name, wrapper)
    (scratch / 'owned_guard.py').symlink_to(GUARD / 'owned_guard.py')
    active_inputs = sorted(contract['source_inputs'] | contract['outputs'])
    active_inputs.extend(str(p) for p in (Path(__file__).resolve(), BASELINE))
    actions = {}
    for stage in STAGES:
        actions[stage] = dict(revision=revision, kind='diagnostic', expected_exit=0,
                              software_provider=False, limits=baseline['replay']['limits'][stage],
                              argv=contract['prospective_emit_commands'][stage.removeprefix('emit-')],
                              frozen_inputs=active_inputs, scope='Pinned ordinary SIMD16 IR replay')
    write_json(scratch / 'actions.json', actions)
    for stage in STAGES:
        check_space()
        print(f'Running {stage}', flush=True)
        result = subprocess.run([sys.executable, '-B', str(wrapper), stage])
        owner_path = scratch / stage / 'owner.json'
        if owner_path.exists():
            owner = read_json(owner_path)
            results[stage] = {key: owner[key] for key in (
                'passed', 'source', 'wall_seconds', 'child_exit', 'stop_reason',
                'peak_sampled_rss_kib', 'whole_inputs_unchanged') if key in owner}
        result.check_returncode()
        product = ir_inventory(Path(actions[stage]['argv'][-1]))
        results[stage]['product'] = product
        if product != baseline['results'][stage]['product']:
            raise RuntimeError(f'Replayed product differed from the pinned {stage} result')


def main():
    receipt, = (Path(arg).resolve() for arg in sys.argv[1:])
    if receipt.exists() or not receipt.is_relative_to(R / 'attestations'):
        raise RuntimeError('Specify a new manifest path under attestations/')
    if Path.cwd() != R or subprocess.check_output(['git', 'status', '--porcelain']).strip():
        raise RuntimeError('Run from the clean Sirius checkout')
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    baseline = read_json(BASELINE)
    inputs = baseline['replay']['source_inputs']
    preset = R / 'bin/linux-gcc'
    bin_created = not preset.parent.exists()
    free = check_space()
    preset.mkdir(parents=True)
    results = {}
    try:
        with tempfile.TemporaryDirectory(prefix='simd16-raw-', dir=R / 'out') as directory:
            scratch = Path(directory)
            recover_sources(inputs)
            acquire_ispc(scratch, inputs)
            contract = generated_sources(inputs, baseline['restored_generated_sources'])
            emit(scratch, revision, baseline, contract, results)
    finally:
        shutil.rmtree(preset)
        if bin_created and not any(preset.parent.iterdir()):
            preset.parent.rmdir()
        if results:
            receipt.parent.mkdir(parents=True, exist_ok=True)
            write_json(receipt, {'scope': 'Ordinary compiler IR replay', 'source': revision,
                                 'driver': pin(Path(__file__)), 'root_free_before_bytes': free,
                                 'baseline': {'path': str(BASELINE.relative_to(R)), **pin(BASELINE)},
                                 'results': results})


if __name__ == '__main__':
    main()
