"""Audit promoted storage and produce the private SIMD16 relocatable object."""

from collections import Counter
from pathlib import Path
import os
import re
import resource
import shlex
import shutil
import subprocess
import sys

import math_source
from resume_raw import BASELINE, GUARD, R, check_space, pin, read_json, write_json
from run_context import PROFILES, tools_at
from run_value_abi import expect_pin, product


SOURCE = Path(__file__).with_name('native_storage.cpp').resolve()
OLD_LOCAL = R / 'out/native-cpu-strict-ir-linking/audit_closed_endpoint_residual.cpp'
OLD_ELF = R / 'out/native-cpu-canonical-closed-remaining/native_check.py'
ELF_PROOF = R / 'out/native-cpu-field-compatibility/elf_proof.py'
ADAPTER_PROOF = Path(__file__).with_name('native_math_check.py').resolve()
MATH_BC = {'bytes': 11580, 'sha256': '60a403f326cdf9c8ada861e7bf2b6a93226b27721909823e4c0ff7ccf8fed2ed'}
MATH_OBJECT = {'bytes': 7400, 'sha256': '919bd088a551000519a1f358d4886db7a08cea801e0319fb47e216d0203b8b3c'}


def derive_local(path):
    if pin(OLD_LOCAL)['sha256'] != '941c668609446d963e8e60381c4b0fa1343fec75d4a8dc9a460ff9598e7d02f4':
        raise RuntimeError('The reviewed Endpoint residual proof changed')
    original = OLD_LOCAL.read_text()
    local = original[original.index('struct LocalCounts '):original.index('struct Audit:')]
    # Scale only the SIMD column bounds and the complete five-column root.
    for before, after in (('O%32', 'O%ColumnBytes'), ('Align(32)', 'Align(ColumnBytes)'),
                          ('uint64_t(O)+32', 'uint64_t(O)+ColumnBytes'), ('Size==32', 'Size==ColumnBytes'),
                          ('Size!=160', 'Size!=5*ColumnBytes'), ('Size==160', 'Size==5*ColumnBytes'),
                          ('!=160', '!=5*ColumnBytes'), ('getNumElements()!=8', 'getNumElements()!=Lanes'),
                          ('sole160B zero', 'sole320B zero')):
        if before not in local:
            raise RuntimeError(f'Endpoint residual source boundary changed: {before}')
        local = local.replace(before, after)
    with path.open('x') as stream:
        stream.write(local)


def native_reader():
    if pin(OLD_ELF)['sha256'] != '99852530e00968fc81f577774e224791fcba4c197e8b8a436ac41914f0c5303e':
        raise RuntimeError('The reviewed ELF/FDE reader changed')
    source = OLD_ELF.read_text()
    if source.count('sirius_local_value_') != 1:
        raise RuntimeError('The named relocation boundary changed')
    source = source.replace('sirius_local_value_', 'sirius_wide_value_')
    namespace = {'__file__': str(OLD_ELF), '__name__': 'retained_wide_native_reader'}
    exec(compile(source, str(OLD_ELF), 'exec'), namespace)
    return namespace['check']


def linked_symbols(path):
    undefined = subprocess.check_output(['/usr/bin/nm', '-u', str(path)], text=True).splitlines()
    external = sorted(line.split()[-1] for line in undefined)
    if set(external) - {'memset', 'sirius_simd_fmaxf'}:
        raise RuntimeError(f'Unexpected external dependency in the partial link: {external}')
    definitions = {line.split()[-1] for line in subprocess.check_output(
        ['/usr/bin/nm', '-g', '--defined-only', str(path)], text=True).splitlines()}
    expected = {f'sirius_wide_strict_{entry}_packet' for entry, _, _ in PROFILES.values()}
    expected.update('sirius_wide_value_' + name + '16' for name in math_source.RECIPES)
    expected.update('sirius_local_value_' + name + '8' for name in math_source.RECIPES)
    if definitions != expected:
        raise RuntimeError(f'Private linked symbol set differs: {sorted(definitions ^ expected)}')
    return {'defined_symbols': sorted(definitions), 'external_symbols': external}


def linked_relocations(path, results):
    expected = Counter()
    for stage in PROFILES:
        expected.update(results['emit-' + stage]['native']['core_relocations'])
    for adapter in results['adapter-object']['adapters'].values():
        expected.update(adapter['core_relocations'])
    relocs = subprocess.check_output(['/usr/bin/readelf', '-Wr', str(path)], text=True)
    actual = Counter(re.findall(
        r'R_X86_64_PLT32\s+[0-9a-f]+\s+(sirius_(?:wide_value|local_value)_\w+)\s+- 4(?:\s|$)',
        relocs))
    symbols = subprocess.check_output(['/usr/bin/readelf', '-Ws', str(path)], text=True)
    hidden = set(re.findall(
        r'^\s*\d+:\s+[0-9a-f]+\s+(?:\d+|0x[0-9a-f]+)\s+FUNC\s+GLOBAL\s+HIDDEN\s+\d+\s+(\w+)$',
        symbols, re.M))
    if actual != expected or hidden != set(results['link']['symbol_audit']['defined_symbols']):
        raise RuntimeError('Partial-link relocations or hidden function bindings differ from the verified objects')
    return {'tool': '/usr/bin/readelf', 'read_flags': ['-Wr', '-Ws'],
            'hidden_global_function_symbols': len(hidden), 'math_plt32_minus4_relocations': dict(sorted(actual.items())),
            'matches_standalone_worker_and_adapter_relocations': True}


def main():
    scratch, context, receipt = (Path(arg).resolve() for arg in sys.argv[1:])
    if not scratch.is_relative_to(R / 'out') or receipt.exists() or not (
        receipt.is_relative_to(scratch) or receipt.is_relative_to(R / 'attestations')
    ):
        raise RuntimeError('Use a task workspace under out/ and a new manifest within the task or attestations/')
    if Path.cwd() != R or subprocess.check_output(['git', 'status', '--porcelain']).strip():
        raise RuntimeError('Run from the clean Sirius checkout')
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    check_space()
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    tools, _ = tools_at(scratch)
    sealed = read_json(context)['results']
    inputs = {stage: scratch / (stage + '.clean.ll') for stage in PROFILES}
    for stage, path in inputs.items():
        record = sealed['cleanup-' + stage]
        if not record['passed'] or record['ir']['owned_context_allocas'] or record['ir']['packing_allocas']:
            raise RuntimeError(f'{stage} did not pass context storage promotion')
        expect_pin(path, record['ir']['output'])
    read_native = native_reader()
    owned = scratch / 'native'
    owned.mkdir()
    include = owned / 'endpoint_local.inc'
    derive_local(include)
    math_ir = R / 'bin/linux-gcc/native-cpu-wide-private/wide-math.ll'
    math_bc, math_o = math_ir.with_suffix('.bc'), math_ir.with_suffix('.o')
    os.link(scratch / 'split-math.ll', math_ir)
    core = math_source.MATH / 'math.o'
    size, hash_value = math_source.PINS['math.o']
    expect_pin(core, {'bytes': size, 'sha256': hash_value})

    flags = shlex.split(subprocess.check_output([str(tools['llvm-config']), '--cxxflags'], text=True))
    plugin = owned / 'sirius-wide-native-storage.so'
    headers = sorted(str(p) for directory in ('sdk/include', 'sdk/lib/clang/23/include')
                     for p in (scratch / directory).rglob('*') if p.is_file())
    commands = {'build': ([str(tools['clang++']), '--driver-mode=g++', *flags, '-O2', '-fPIC', '-shared',
                           '-I' + str(owned), str(SOURCE), '-o', str(plugin)], headers)}
    for stage, path in inputs.items():
        commands['audit-' + stage] = ([str(tools['opt']), '-load-pass-plugin=' + str(plugin),
            '-passes=sirius-wide-native-storage,verify', '-verify-each', '-disable-output', str(path)],
            [str(plugin), str(path)])
    commands['adapter-ir'] = ([str(tools['opt']), '-passes=verify', '-verify-each', str(math_ir), '-o', str(math_bc)], [str(math_ir)])
    def lower(path, output, cpu):
        return [str(tools['llc']), '-O2', '-filetype=obj', '-relocation-model=pic', '-mcpu=' + cpu,
                str(path), '-o', str(output)]

    commands['adapter-object'] = (lower(math_bc, math_o, 'x86-64-v4'), [str(math_bc)])
    objects = [str(math_o), str(core)]
    for stage, path in inputs.items():
        output = owned / (stage + '.o')
        commands['emit-' + stage] = (lower(path, output, 'skylake-avx512'), [str(path)])
        objects.append(str(output))
    commands['link'] = (['/usr/bin/ld', '-r', *objects, '-o', str(owned / 'simd16-component.o')], objects)

    shutil.copyfile(GUARD / 'run_owned.py', owned / 'run_owned.py')
    (owned / 'owned_guard.py').symlink_to(GUARD / 'owned_guard.py')
    common = [str(p) for p in (context, SOURCE, OLD_LOCAL, OLD_ELF, ELF_PROOF, ADAPTER_PROOF, include,
                               Path(__file__).resolve(), *tools.values(), Path('/usr/bin/ld'), Path('/usr/bin/nm'),
                               Path('/usr/bin/objdump'), Path('/usr/bin/readelf'))]
    limits = read_json(BASELINE)['replay']['limits']['emit-transport']
    write_json(owned / 'actions.json', {name: dict(revision=revision, kind='diagnostic', expected_exit=0,
        software_provider=False, limits=limits, argv=argv, frozen_inputs=common + dependencies,
        build_environment={'TMPDIR': str(scratch)}, scope='SIMD16 private native storage and component gate')
        for name, (argv, dependencies) in commands.items()})

    results = {}
    try:
        for name, (argv, _) in commands.items():
            check_space()
            status = subprocess.run([sys.executable, '-B', str(owned / 'run_owned.py'), name])
            owner_file = owned / name / 'owner.json'
            log = owned / name / 'stderr.log'
            if owner_file.exists():
                owner = read_json(owner_file)
                results[name] = {key: owner[key] for key in ('passed', 'source', 'wall_seconds', 'child_exit',
                    'stop_reason', 'whole_inputs_unchanged', 'peak_sampled_rss_kib') if key in owner}
                if not owner['passed']:
                    results[name]['diagnostic'] = log.read_text()[:1500]
            status.check_returncode()
            result = results[name]
            output = Path(argv[-1])
            if name.startswith('audit-'):
                summary = re.search(r'^SIMD16 native storage: \w+ staging_bytes=(\d+) residual_allocas=(\d+) residual_bytes=(\d+)$',
                                    log.read_text(), re.M)
                if not summary: raise RuntimeError(f'Missing read-only storage proof for {name}')
                result['storage'] = dict(zip(('staging_bytes', 'residual_allocas', 'residual_bytes'), map(int, summary.groups())))
                continue
            result['output'] = product(output)
            if name == 'adapter-ir':
                expect_pin(output, MATH_BC)
            elif name == 'adapter-object':
                expect_pin(output, MATH_OBJECT)
                from native_math_check import all_adapters
                names = ['sirius_wide_value_' + recipe + '16' for recipe in math_source.RECIPES]
                checked = all_adapters(output, names)
                for recipe, symbol in zip(math_source.RECIPES, names):
                    if checked[symbol]['core_relocations'] != {'sirius_local_value_' + recipe + '8': 2}:
                        raise RuntimeError(f'Adapter core call join differs: {recipe}')
                result['adapters'] = checked
            elif name.startswith('emit-'):
                stage = name.removeprefix('emit-')
                entry = sealed['cleanup-' + stage]['ir']['entry']
                checked = read_native(output, entry)
                if set(checked['core_relocations']) != set(sealed['cleanup-' + stage]['ir']['value_calls']):
                    raise RuntimeError(f'{stage} native value target set differs')
                result['native'] = checked
            elif name == 'link':
                result['symbol_audit'] = linked_symbols(output)
                result['relocation_audit'] = linked_relocations(output, results)
    finally:
        if results:
            receipt.parent.mkdir(parents=True, exist_ok=True)
            write_json(receipt, {'scope': 'Compiler, standalone ELF and relocatable component only', 'source': revision,
                'driver': product(Path(__file__).resolve()), 'storage_pass': product(SOURCE),
                'context': product(context), 'endpoint_proof_source': product(include), 'results': results})


if __name__ == '__main__':
    main()
