"""Lower and verify the initialized context of the three sealed value modules."""

from collections import Counter
from pathlib import Path
import hashlib
import re
import resource
import shlex
import shutil
import subprocess
import sys

from recover_sdk_ranges import descriptor, pinned_tools
from resume_raw import BASELINE, GUARD, R, check_space, read_json, write_json
from run_value_abi import product
from transport_context import BASE, derive


SEALED = R / 'attestations/qualification/432d65a9/native-cpu-simd16-value/MANIFEST.json'
SOURCE = Path(__file__).with_name('context_lowering.cpp').resolve()
PROFILES = {'transport': ('transport', 4699, 184384), 'endpoint': ('endpoint', 3829, 235584),
            'dopri_phase': ('dp', 4149, 241664)}
CLEANUP = 'function(sroa,mem2reg,instcombine,simplifycfg,early-cse<memssa>,dse),globaldce,verify'


def tools_at(scratch):
    sdk, catalog, _ = descriptor()
    return pinned_tools(scratch / 'sdk', sdk, ('clang++', 'llvm-config', 'opt', 'llc')), catalog


def value_identity(path, sealed_output, sealed_strict):
    previous = f"; ModuleID = '{R / sealed_strict['path']}'\n".encode()
    current = f"; ModuleID = '{path.with_name(path.name.replace('.value.ll', '.strict.ll'))}'\n".encode()
    with path.open('rb') as stream:
        header = stream.readline()
        if header not in (previous, current) or path.stat().st_size - len(header) + len(previous) != sealed_output['bytes']:
            raise RuntimeError(f'Value replay header or size differs: {path}')
        digest = hashlib.sha256(previous)
        while data := stream.read(1024 * 1024): digest.update(data)
    if digest.hexdigest() != sealed_output['sha256']:
        raise RuntimeError(f'Value replay differs beyond the ModuleID filename: {path}')
    return {'input': product(path), 'equals_sealed_with_moduleid_restored': True}


def summarize(path, stage, cleaned):
    text = path.read_text()
    entry, words, _ = PROFILES[stage]
    symbol = f'sirius_wide_strict_{entry}_packet'
    headers = re.findall(r'^define .*$', text, re.M)
    if len(headers) != 1 or f'@{symbol}(' not in headers[0] or 'hidden' not in headers[0].split():
        raise RuntimeError(f'{stage} hidden packet entry differs')
    attrs = dict(re.findall(r'^attributes #(\d+) = (.*)$', text, re.M))
    attr = attrs[re.search(r'#(\d+)', headers[0])[1]]
    if 'strictfp' not in attr or 'uwtable' not in attr or 'uwtable(sync)' in attr:
        raise RuntimeError(f'{stage} strict FP or asynchronous unwind attribute differs')
    if re.search(r'= (?:fadd|fsub|fmul|fdiv|frem|fcmp|sitofp|uitofp|fptosi|fptoui|fpext|fptrunc)\b', text):
        raise RuntimeError(f'{stage} ordinary FP arithmetic was introduced')
    constrained = re.findall(r'\bcall\s+[^\n]*?@(llvm\.experimental\.constrained\.[\w.$-]+)\([^\n]*', text)
    for operation in re.findall(r'\bcall\s+[^\n]*?@llvm\.experimental\.constrained\.[^\n]+', text):
        if 'metadata !"fpexcept.strict"' not in operation or (
            '@llvm.experimental.constrained.fcmp' not in operation and 'metadata !"round.dynamic"' not in operation
        ):
            raise RuntimeError(f'{stage} constrained FP mode changed')

    allocas = re.findall(r'^  %[^=\n]+ = alloca .+$', text, re.M)
    staging = [line for line in allocas if line.startswith('  %staged_output =')]
    if staging != [f'  %staged_output = alloca [{words * 16} x i32], align 64']:
        raise RuntimeError(f'{stage} transactional staging differs')
    context = [line for line in allocas if 'alloca %v16_varying_KernelContext_0' in line]
    banks = [line for line in allocas if re.search(r'alloca \[(?:240|48) x i32\]', line) or
             line.startswith(('  %SiriusBank.', '  %SiriusPackedBank.', '  %SiriusPackedTags.'))]
    calls = Counter(re.findall(r'\bcall x86_regcallcc [^\n]*?@(sirius_wide_value_\w+)\(', text))
    if not calls or 'sirius_wide_strict_pointer_' in text:
        raise RuntimeError(f'{stage} lost the strict value ABI')
    return {'output': product(path), 'entry': symbol, 'staging_bytes': words * 16 * 4,
            'owned_context_allocas': len(context), 'packing_allocas': len(banks), 'value_calls': dict(calls),
            'constrained_fp_calls': dict(Counter(constrained)),
            'alloca_count': len(allocas), 'residual_allocas': [line.strip() for line in allocas if line not in staging]
            if cleaned and len(allocas) < 32 else []}


def main():
    scratch, receipt = (Path(arg).resolve() for arg in sys.argv[1:])
    if not scratch.is_relative_to(R / 'out') or receipt.exists() or not (
        receipt.is_relative_to(scratch) or receipt.is_relative_to(R / 'attestations')
    ):
        raise RuntimeError('Use a task workspace under out/ and a new manifest within the task or attestations/')
    if Path.cwd() != R or subprocess.check_output(['git', 'status', '--porcelain']).strip():
        raise RuntimeError('Run from the clean Sirius checkout')
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    check_space()
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    tools, catalog = tools_at(scratch)
    sealed = read_json(SEALED)['results']
    raw = {stage: scratch / (stage + '.value.ll') for stage in PROFILES}
    admitted = {stage: value_identity(path, sealed['value-' + stage]['output'], sealed['scope-' + stage]['output'])
                for stage, path in raw.items()}
    owned = scratch / 'context'
    owned.mkdir()
    include = owned / 'transport_context.inc'
    derive(include)
    support = scratch / 'sha256-only.a'
    flags = shlex.split(subprocess.check_output([str(tools['llvm-config']), '--cxxflags'], text=True))
    plugin = owned / 'sirius-wide-context.so'
    headers = sorted(str(p) for directory in ('sdk/include', 'sdk/lib/clang/23/include')
                     for p in (scratch / directory).rglob('*') if p.is_file())
    commands = {'build': ([str(tools['clang++']), '--driver-mode=g++', *flags, '-O2', '-fPIC',
                          '-shared', '-Wno-deprecated-declarations', '-I' + str(owned), str(SOURCE), str(support),
                          '-Wl,--exclude-libs,ALL', '-o', str(plugin)], [str(support), *headers])}
    for stage, path in raw.items():
        normal = scratch / (stage + '.context.ll')
        commands['context-' + stage] = ([str(tools['opt']), '-load-pass-plugin=' + str(plugin),
            '-passes=sirius-wide-context,verify', '-verify-each', '-S', str(path), '-o', str(normal)],
            [str(plugin), str(path)])
    for stage in raw:
        normal = scratch / (stage + '.context.ll')
        commands['cleanup-' + stage] = ([str(tools['opt']), '-passes=' + CLEANUP, '-mcpu=skylake-avx512',
            '-verify-each', '-S', str(normal), '-o', str(scratch / (stage + '.clean.ll'))], [str(normal)])
    shutil.copyfile(GUARD / 'run_owned.py', owned / 'run_owned.py')
    (owned / 'owned_guard.py').symlink_to(GUARD / 'owned_guard.py')
    common = [str(path) for path in (SEALED, SOURCE, BASE, include, Path(__file__).resolve(), *tools.values())]
    limits = read_json(BASELINE)['replay']['limits']['emit-transport']
    write_json(owned / 'actions.json', {name: dict(revision=revision, kind='diagnostic', expected_exit=0,
        software_provider=False, limits=limits, argv=argv, frozen_inputs=common + keys,
        build_environment={'TMPDIR': str(scratch)}, scope='Pinned SIMD16 context compiler gate')
        for name, (argv, keys) in commands.items()})
    results = {}
    try:
        for name, (argv, _) in commands.items():
            check_space()
            status = subprocess.run([sys.executable, '-B', str(owned / 'run_owned.py'), name])
            owner_file = owned / name / 'owner.json'
            if owner_file.exists():
                owner = read_json(owner_file)
                results[name] = {key: owner[key] for key in ('passed', 'source', 'wall_seconds', 'child_exit',
                    'stop_reason', 'whole_inputs_unchanged', 'peak_sampled_rss_kib') if key in owner}
                if not owner['passed']:
                    results[name]['diagnostic'] = (owned / name / 'stderr.log').read_text()[:1500]
            status.check_returncode()
            output = Path(argv[-1])
            if name == 'build':
                results[name]['output'] = product(output)
                continue
            phase, stage = name.split('-', 1)
            cleaned = phase == 'cleanup'
            detail = summarize(output, stage, cleaned)
            results[name]['ir'] = detail
            if cleaned and (detail['owned_context_allocas'] or detail['packing_allocas']):
                raise RuntimeError(f'{stage} main context or packing bank survived cleanup')
            if not cleaned:
                log = (owned / name / 'stderr.log').read_text()
                entry, words, size = PROFILES[stage]
                match = re.search(r'^SIMD16 context: ' + entry +
                    r' context_bytes=(\d+) staging_words=(\d+) bounded_accesses=(\d+) normalized=(\d+) strict_operations=(\d+)$',
                    log, re.M)
                if not match or tuple(map(int, match.groups()[:2])) != (size, words * 16):
                    raise RuntimeError(f'{stage} context proof inventory disagrees with the layout')
                detail.update(context_bytes=size, bounded_private_accesses=int(match[3]),
                              normalized_private_accesses=int(match[4]), retained_strict_operations=int(match[5]))
    finally:
        if results:
            receipt.parent.mkdir(parents=True, exist_ok=True)
            write_json(receipt, {'scope': 'SIMD16 context normalization and compiler storage promotion',
                                 'source': revision, 'driver': product(Path(__file__).resolve()),
                                 'pass': product(SOURCE), 'transport_proof': product(include),
                                 'sealed_value': product(SEALED), 'sdk_archive': catalog['archive'],
                                 'value_replays': admitted, 'results': results})


if __name__ == '__main__':
    main()
