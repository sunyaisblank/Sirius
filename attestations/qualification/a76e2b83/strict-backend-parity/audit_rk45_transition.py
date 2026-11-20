"""Isolate the historical RK45 rejection change through its two real tracer objects."""

import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess


ROOT = Path.cwd()
TASK = ROOT / 'out/strict-backend-science'
WORK = TASK / 'historical-rk45'
SOURCE = WORK / 'src'
PROBE = Path(__file__).with_name('historical-rk45-probe.cpp').resolve()
COMPILER = '/usr/bin/g++-14'
FLAGS = ['-O3', '-DNDEBUG', '-std=c++2c', '-Wall', '-Wextra', '-Wpedantic',
         '-Wno-missing-field-initializers', '-Werror', '-ffp-contract=off']
AFTER = 'e736ffe0c318dcd637f1453beb3523ad14a15047'
TRACER = 'sirius/backend/cpu/geodesic_tracer.cpp'
SHARED = ('sirius/base/contracts.cpp', 'sirius/core/tensor.cpp', 'sirius/core/geodesic_integrator.cpp')


def pin(path):
    with path.open('rb') as source:
        sha = hashlib.file_digest(source, 'sha256').hexdigest()
    return {'path': str(path.relative_to(ROOT)), 'bytes': path.stat().st_size, 'sha256': sha}


def run():
    previous = subprocess.check_output(['git', 'rev-parse', AFTER + '^'], text=True).strip()
    WORK.mkdir(parents=True, exist_ok=True)
    if not SOURCE.exists():
        export = subprocess.check_output(['git', 'archive', '--format=tar', AFTER, 'src/sirius'])
        subprocess.run(['/usr/bin/tar', '-x', '-C', str(WORK)], input=export, check=True)
    original_tracer = WORK / 'geodesic_tracer_before.cpp'
    if not original_tracer.exists():
        original_tracer.write_bytes(subprocess.check_output(['git', 'cat-file', 'blob',
                                                               f'{previous}:src/{TRACER}']))
    compile_commands = []

    def compile_object(source, output, extra=()):
        argv = [COMPILER, *FLAGS, '-I', str(SOURCE), *extra, '-c', str(source), '-o', str(output)]
        subprocess.run(argv, check=True, timeout=180)
        compile_commands.append(argv)
        return output

    compare_paths = subprocess.check_output(['git', 'diff', '--name-only', previous, AFTER, '--',
                                             'src/sirius/backend/cpu', 'src/sirius/core',
                                             'src/sirius/base'], text=True).splitlines()
    if set(compare_paths) != {'src/' + TRACER, 'src/sirius/core/metrics/morris_thorne_family.h',
                            'src/sirius/core/metrics/registry.h'}:
        raise RuntimeError('Unexpected historical numerical source changes')
    tracer_sources = {'before': original_tracer, 'after': SOURCE / TRACER}
    for revision, path in zip((previous, AFTER), tracer_sources.values(), strict=True):
        if path.read_bytes() != subprocess.check_output(['git', 'cat-file', 'blob',
                                                        f'{revision}:src/{TRACER}']):
            raise RuntimeError(f'The recovered tracer differs from {revision}')

    shared = [compile_object(SOURCE / name, WORK / (Path(name).stem + '.o')) for name in SHARED]
    objects = {variant: compile_object(source, WORK / f'{variant}.o')
               for variant, source in tracer_sources.items()}
    undefined = subprocess.check_output(['/usr/bin/nm', '-u', str(objects['before'])], text=True)
    symbol, = (line.split()[-1] for line in undefined.splitlines() if 'IntegrateStepRk45' in line)
    probe = compile_object(PROBE, WORK / 'probe.o', [f'-DSIRIUS_RK45_LINK_SYMBOL="{symbol}"'])
    rows = {}
    inputs = {}
    for variant, tracer in objects.items():
        executable = WORK / variant
        argv = [COMPILER, f'-Wl,--wrap={symbol}', *map(str, (probe, tracer, *shared)), '-o', str(executable)]
        subprocess.run(argv, check=True, timeout=180)
        captured = subprocess.check_output([executable], text=True, timeout=60)
        inputs[variant] = {'tracer': pin(tracer_sources[variant]), 'executable': pin(executable),
                           'link_command': argv, 'trace_transcript_sha256':
                               hashlib.sha256(captured.encode()).hexdigest()}
        rows[variant] = {(row['profile'], row['x'], row['y']): row
                         for row in csv.DictReader(io.StringIO(captured))}

    if rows['before'].keys() != rows['after'].keys() or not rows['after']:
        raise RuntimeError('The two traces did not cover the same witness rays')
    witnesses = {'default_step': [], 'rejection_stress': []}
    control = {'captured': 0, 'escaped': 0}
    for key, before in rows['before'].items():
        after = rows['after'][key]
        if before['impact'] != after['impact'] or before['analytic'] != after['analytic']:
            raise RuntimeError(f'Initial phase differs on ray {key}')
        if before['changed_phase'] != '0' or after['changed_phase'] != '0':
            raise RuntimeError(f'A retry changed the committed phase on ray {key}')
        if after['outcome'] != after['analytic'] or after['numerical_failure'] != '0':
            raise RuntimeError(f'The revised tracer disagrees with the analytic classification on ray {key}')
        control[after['outcome']] += 1
        if (before['outcome'] == 'max_steps' and before['rejections'] == '1' and
                before['numerical_failure'] == '0' and int(before['attempts']) < 5000 and
                int(after['attempts']) > int(before['attempts'])):
            witnesses[key[0]].append({'before': before, 'after': after})

    if not witnesses['rejection_stress'] or not all(control.values()):
        raise RuntimeError('The diagnostic did not witness a recoverable rejection and both physical sides')
    data = {'scope': 'NEW_HISTORICAL_RAY_DIAGNOSTIC', 'before_revision': previous,
            'after_revision': AFTER, 'full_image_reference_or_July_tape_claimed': False,
            'exported_source_tree': subprocess.check_output(
                ['git', 'rev-parse', f'{AFTER}:src/sirius'], text=True).strip(),
            'isolated_changed_translation_unit': 'src/' + TRACER,
            'compiler_version': subprocess.check_output([COMPILER, '-dumpfullversion'], text=True).strip(),
            'analytic_oracle': 'The Schwarzschild null angular-momentum barrier (1-2/r)L^2/r^2 peaks at r=3 for M=1. '
                               'Incoming rays from r=50 capture below L/E=sqrt(27), escape above it.',
            'minimum_impact_clearance': 0.2, 'maximum_trace_attempts': 5000,
            'probe': pin(PROBE), 'shared_objects': [pin(path) for path in shared],
            'compile_commands': compile_commands, 'variants': inputs,
            'total_rays_per_variant': len(rows['after']), 'revised_outcomes': control,
            'recoverable_rejection_witnesses': witnesses,
            'all_rows': {variant: list(group.values()) for variant, group in rows.items()}}
    (TASK / 'historical-rk45-result.json').write_text(json.dumps(data, sort_keys=True, indent=2) + '\n')
    print({name: len(group) for name, group in witnesses.items()}, control)


if __name__ == '__main__':
    run()
