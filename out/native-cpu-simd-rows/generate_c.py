"""Seal, then generate the three original row-owner shaders as C once."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
BASIS = ROOT / 'out/native-cpu-slang-bridge'
OUTPUT = ROOT / 'bin/linux-gcc/native-cpu-simd-rows/c-source'
PREPARED = WORK / 'prepared-c-generation.json'
RESULT = WORK / 'c-generation-result.json'
PINS = {
    'out/native-cpu-row-ownership/compiler-result.json': '05b5a01dcaddeba1deacca9ed73af0ee64c1d20cc512d2ad90948aa987a50449',
    'out/native-cpu-slang-bridge/batch-generation.json': 'de6d68848c352c3b9262e01822730d25459b725d42f57e4b1769ba93d2cc2428',
    'out/native-cpu-slang-bridge/generation-result.json': '57e149a7ee8924e49f64ea39383502c6ce017ab84e5a64c3bfbf72d8481c849b',
    'out/native-cpu-slang-bridge/primitive-lowering.patch': 'c7eabc7957b1b26dfb383ebd1781f8cb5edc3414e6bed968de9c119cf6142d8c',
    'out/native-cpu-slang-bridge/generate_cpp.py': '7fb8444fc0e8af7585df5d4a895e18806ea3ce3b5ccb65ede9247eb2f6f92d44',
    'out/native-cpu-slang-bridge/batch_generate.py': '7dcb6ad2c97b6edb1e717169e9e65e8536b016272d4448f7a8d32fea307f4900',
    'out/native-cpu-slang-bridge/prepare_batch.py': '5aa607cefd69d77bb033146ceeb539dffc7152fa5189f7744cc694753dbd7d68',
    'out/native-cpu-slang-bridge/build/actions.json': 'a31728c3cd50b9651864c7a15fba65382072ed26f51bb07f76bcd39b69f6d67f',
}
PROGRAMS = {
    'transport': ('9c581d3a98843f9cd831a7814e8b7ad77a0f40ab2c8cb29301df67cc6bad1ffc', 15967),
    'endpoint': ('b6f8ba278b7e9df65860143d9711972603b4f8a1850b2ade138bf17f625150d5', 21143),
    'dopri_phase': ('3df9a067a4e65fd58959fd84076bb6c743e33fb15f26aa898d699029c6cd64fe', 10689),
}


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def clean_revision():
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    return subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()


def prepare():
    revision = clean_revision()
    assert not any(p.exists() for p in (PREPARED, RESULT, OUTPUT, WORK / 'generation/actions.json'))
    for name, digest in PINS.items():
        assert pin(ROOT / name)['sha256'] == digest, name
    authority = read(ROOT / 'out/native-cpu-row-ownership/compiler-result.json')
    batch = read(BASIS / 'batch-generation.json')
    dp = read(BASIS / 'generation-result.json')
    old_actions = read(BASIS / 'build/actions.json')
    names = set(PINS)
    latest = read(BASIS / 'build/batch-source/inputs-after.json')
    for action_name in ('batch-source', 'dopri-cpp-native-primitives'):
        old = old_actions[action_name]
        # Existing owned receipts bind the actual original compiler and imports.
        saved = read(BASIS / 'build' / action_name / 'inputs-after.json')
        for name in old['frozen_inputs']:
            path = ROOT / name
            current = pin(path)
            # The later three-stage authority supersedes the initial DP-only
            # source contract; retain the older receipt for unique DP inputs.
            expected = latest.get(name, saved[name])
            assert current == {k: expected[k] for k in ('bytes', 'sha256')}, name
            names.add(name)
    stages = []
    for kind, slug in (('Transport', 'transport'), ('Endpoint', 'endpoint'), ('DopriPhase', 'dopri_phase')):
        original = list((dp if kind == 'DopriPhase' else batch['modules'][kind])['actual_CPP_argv'])
        assert original[0] == authority['tools']['compiler']
        assert original[1] == str(BASIS / 'shaders' / ('retained_' + slug + '.slang'))
        assert original.count('-DSIRIUS_RETAINED_CPU_CPP=1') == 1
        assert original[original.index('-target') + 1] == 'cpp'
        assert original[original.index('-fp-mode') + 1] == 'precise'
        assert '-profile' not in original
        argv = list(original)
        source = OUTPUT / (slug + '.c')
        reflection = OUTPUT / (slug + '.reflection.json')
        substitutions = {original.index('-target') + 1: 'c', original.index('-o') + 1: str(source),
                         original.index('-reflection-json') + 1: str(reflection)}
        for index, value in substitutions.items():
            argv[index] = value
        assert {i for i, (a, b) in enumerate(zip(original, argv)) if a != b} == set(substitutions)
        module = authority['modules'][kind]
        assert pin(ROOT / module['module']) == {k: module[k] for k in ('bytes', 'sha256')}
        names.add(module['module'])
        table = ROOT / 'bin/linux-gcc/native-cpu-slang-bridge/typed-batch' / (slug + '.program.hpp')
        assert pin(table)['sha256'] == PROGRAMS[slug][0]
        names.add(str(table.relative_to(ROOT)))
        stages.append({'stage': kind, 'original_CPP_argv': original, 'argv': argv,
                       'source': str(source.relative_to(ROOT)), 'reflection': str(reflection.relative_to(ROOT)),
                       'original_SPIRV': module, 'canonical_program': {'path': str(table.relative_to(ROOT)),
                       **pin(table), 'words': PROGRAMS[slug][1]}})
    names |= {str(p.relative_to(ROOT)) for p in (BASIS / 'shaders').iterdir() if p.is_file()}
    names |= {str((WORK / name).relative_to(ROOT)) for name in
              ('generate_c.py', 'generation/run_owned.py', 'generation/owned_guard.py')}
    inputs = {name: pin(ROOT / name) for name in sorted(names)}
    prepared = {'revision': revision, 'inputs': inputs, 'stages': stages,
                'scope': 'Original three row-owner shaders; existing four primitive lowerings unchanged. Only CPP-to-C target and output paths differ. No SPIRV, host, native object or scientific workload generation.'}
    write(PREPARED, prepared)
    action = {**old_actions['batch-source'], 'revision': revision,
              'argv': ['/usr/bin/python3', '-B', str(WORK / 'generate_c.py'), 'run'],
              'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
              'frozen_inputs': sorted(names | {str(PREPARED.relative_to(ROOT))}),
              'scope': prepared['scope']}
    write(WORK / 'generation/actions.json', {'c-source': action})
    print(json.dumps({'prepared_generations': 3, 'generated': False, 'revision': revision}))


def run():
    prepared = read(PREPARED)
    assert clean_revision() == prepared['revision']
    assert not RESULT.exists() and not OUTPUT.exists()
    for name, expected in prepared['inputs'].items():
        assert pin(ROOT / name) == expected, name
    OUTPUT.mkdir(parents=True)
    record = {'revision': prepared['revision'], 'scope': prepared['scope'], 'stages': [],
              'all_generations_succeeded': False, 'shader_execution': False}
    write(RESULT, record)
    for stage in prepared['stages']:
        row = {**stage, 'status': 'started'}
        record['stages'].append(row)
        write(RESULT, record)
        print(json.dumps({'stage': stage['stage'], 'argv': stage['argv']}), flush=True)
        result = subprocess.run(stage['argv'], cwd=ROOT)
        row.update(status='terminal', exit=result.returncode)
        for key in ('source', 'reflection'):
            path = ROOT / stage[key]
            if path.is_file():
                row[key + '_identity'] = pin(path)
        write(RESULT, record)
    record['inputs_unchanged'] = all(pin(ROOT / name) == value for name, value in prepared['inputs'].items())
    record['source_unchanged'] = clean_revision() == prepared['revision']
    record['all_generations_succeeded'] = all(row['exit'] == 0 and row.get('source_identity', {}).get('bytes', 0) > 0
                                              for row in record['stages'])
    write(RESULT, record)
    assert record['inputs_unchanged'] and record['source_unchanged']
    raise SystemExit(0 if record['all_generations_succeeded'] else 1)


if __name__ == '__main__':
    assert len(sys.argv) == 2 and sys.argv[1] in ('prepare', 'run')
    {'prepare': prepare, 'run': run}[sys.argv[1]]()
