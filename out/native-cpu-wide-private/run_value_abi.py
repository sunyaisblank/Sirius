"""Replay the sealed strict inputs and verify their SIMD16 value ABI."""

from collections import Counter
from pathlib import Path
import re
import resource
import shutil
import subprocess
import sys
import zipfile

import math_source
from recover_sdk_ranges import CATALOG, pinned_tools, recover
import run_strict_scope as strict
from resume_raw import BASELINE, GUARD, R, check_space, pin, read_json, restore, write_json


SCOPE = R / 'attestations/qualification/f292c239/native-cpu-simd16-strict/MANIFEST.json'
MATH_CAPSULE = R / 'attestations/qualification/b2cbd4db/native-cpu-private-context-interface-gates'
VALUE_SOURCE = Path(__file__).with_name('value_abi.cpp').resolve()
VALUE_PREFIX = 'sirius_wide_value_'
POINTER_PREFIX = 'sirius_wide_strict_pointer_'


def product(path):
    return {'path': str(path.relative_to(R)), **pin(path)}


def expect_pin(path, record):
    if pin(path) != {key: record[key] for key in ('bytes', 'sha256')}:
        raise RuntimeError(f'Pinned product differs: {path.relative_to(R)}')


def recover_adapter(scratch):
    for bundle_name, names in (('source-and-evidence.zip', ('math.ll',)), ('products.zip', ('math.bc', 'math.o'))):
        with zipfile.ZipFile(MATH_CAPSULE / bundle_name) as bundle:
            for name in names:
                path = math_source.MATH / name
                size, sha = math_source.PINS[name]
                expected = {'bytes': size, 'sha256': sha}
                if path.exists():
                    expect_pin(path, expected)
                    continue
                with bundle.open(str(path.relative_to(R))) as stream:
                    restore(stream, path, expected)
    text, _ = math_source.sources()
    adapter = scratch / 'split-math.ll'
    with adapter.open('x') as stream:
        stream.write(text)
    return adapter


def headers(text, kind):
    # LLVM prints the aggregate and parameter types on one line. The defined
    # adapter has argument names. Stage declarations carry the same types.
    result = {}
    for line in text.splitlines():
        if not line.startswith(kind + ' ') or '@' + VALUE_PREFIX not in line:
            continue
        match = re.fullmatch(r'(?:declare|define) (?:dso_local )?hidden x86_regcallcc (\{ [^}]+ \}) @('
                             + VALUE_PREFIX + r'\w+)\(([^)]*)\).*', line)
        if not match:
            raise RuntimeError(f'Unexpected value ABI declaration: {line}')
        result[match[2]] = (match[1], tuple(param.split(' %', 1)[0] for param in match[3].split(', ')))
    return result


def audit_stage(stage, path, strict_path, adapter_headers, audited):
    before = strict_path.read_text()
    after = path.read_text()
    expected = Counter(re.findall(r'\bcall void @' + POINTER_PREFIX + r'(\w+)\(', before))
    observed = Counter(re.findall(r'\bcall x86_regcallcc [^\n]*?@' + VALUE_PREFIX + r'(\w+)\(', after))
    if observed != expected or sum(expected.values()) != audited['checked_pointer_calls'] or POINTER_PREFIX in after:
        raise RuntimeError(f'{stage} pointer-to-value call migration is incomplete')
    stage_headers = headers(after, 'declare')
    if set(stage_headers) != {VALUE_PREFIX + name for name in expected} or any(
        adapter_headers[name] != value for name, value in stage_headers.items()
    ):
        raise RuntimeError(f'{stage} value declarations differ from the verified split adapter')
    fp_pattern = r'\bcall\s+[^\n]*?@(llvm\.experimental\.constrained\.[\w.$-]+|llvm\.fabs\.\w+|sirius_simd_fmaxf)\('
    if re.findall(fp_pattern, before) != re.findall(fp_pattern, after):
        raise RuntimeError(f'{stage} strict FP call sequence changed')
    if len(re.findall(r'= fneg ', after)) != audited['retained_fneg'] or re.search(
        r'= (?:fadd|fsub|fmul|fdiv|frem|fcmp|sitofp|uitofp|fptosi|fptoui|fpext|fptrunc)\b', after
    ):
        raise RuntimeError(f'{stage} ordinary FP arithmetic changed')
    context = r'^%v16_varying_KernelContext_0 = type .+$'
    if re.search(context, before, re.M)[0] != re.search(context, after, re.M)[0]:
        raise RuntimeError(f'{stage} context layout changed during value binding')
    return {'output': product(path), 'value_calls': dict(observed), 'adapter_declarations_match': True,
            'strict_fp_sequence_preserved': True, 'context_layout_preserved': True}


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
    sealed = read_json(SCOPE)['results']
    limits = read_json(BASELINE)['replay']['limits']['emit-transport']
    raw = {stage: R / sealed[stage]['raw']['path'] for stage in strict.STAGES}
    for stage, path in raw.items():
        expect_pin(path, sealed[stage]['raw'])

    descriptor = read_json(strict.HISTORICAL_SDK)['setup']['official_sdk']
    sdk_directory = scratch / 'sdk'
    if not sdk_directory.exists():
        recover(scratch)
    tools = pinned_tools(sdk_directory, descriptor, strict.TOOLS)
    adapter = recover_adapter(scratch)
    verified_adapter = scratch / 'split-math.checked.ll'
    strict_build = strict.build_action(scratch, tools)
    value_build = [str(VALUE_SOURCE) if arg == str(strict.SOURCE) else arg for arg in strict_build]
    value_build[-1] = str(scratch / 'sirius-wide-value.so')
    owned = scratch / 'value'
    owned.mkdir()
    shutil.copyfile(GUARD / 'run_owned.py', owned / 'run_owned.py')
    (owned / 'owned_guard.py').symlink_to(GUARD / 'owned_guard.py')
    header_paths = sorted(str(path) for directory in ('sdk/include', 'sdk/lib/clang/23/include')
                          for path in (scratch / directory).rglob('*') if path.is_file())
    common = [str(path) for path in (SCOPE, strict.HISTORICAL_SDK, BASELINE, CATALOG, VALUE_SOURCE,
                                     strict.SOURCE, Path(__file__).resolve(), *tools.values())]
    build_inputs = [*header_paths, strict_build[-4]]
    commands = {
        'build-scope': (strict_build, build_inputs),
        'build-value': (value_build, build_inputs),
        'adapter': ([str(tools['opt']), '-passes=verify', '-verify-each', '-S', str(adapter), '-o', str(verified_adapter)],
                    [str(adapter), *[str(math_source.MATH / name) for name in math_source.PINS]]),
    }
    for stage, path in raw.items():
        strict_ir = scratch / (stage + '.strict.ll')
        commands['scope-' + stage] = ([str(tools['opt']), '-load-pass-plugin=' + strict_build[-1],
            '-passes=sirius-wide-strict,verify', '-verify-each', '-S', str(path), '-o', str(strict_ir)],
            [strict_build[-1], str(path)])
        commands['value-' + stage] = ([str(tools['opt']), '-load-pass-plugin=' + value_build[-1],
            '-passes=sirius-wide-value,verify', '-verify-each', '-S', str(strict_ir),
            '-o', str(scratch / (stage + '.value.ll'))], [value_build[-1], str(strict_ir), str(verified_adapter)])

    actions = {name: dict(revision=revision, kind='diagnostic', expected_exit=0, software_provider=False,
                         limits=limits, argv=argv, frozen_inputs=common + inputs,
                         build_environment={'TMPDIR': str(scratch)}, scope='Pinned SIMD16 value ABI compiler gate')
               for name, (argv, inputs) in commands.items()}
    write_json(owned / 'actions.json', actions)
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
                    results[name]['diagnostic'] = (owned / name / 'stderr.log').read_text()[:1200]
            status.check_returncode()
            output = Path(argv[-1])
            if name == 'adapter':
                adapter_headers = headers(output.read_text(), 'define')
                if set(adapter_headers) != {VALUE_PREFIX + recipe + '16' for recipe in math_source.RECIPES}:
                    raise RuntimeError('LLVM adapter entry points differ')
            elif name.startswith('scope-'):
                stage = name.removeprefix('scope-')
                expect_pin(output, sealed[stage]['scope']['output'])
            elif name.startswith('value-'):
                stage = name.removeprefix('value-')
                detail = audit_stage(stage, output, scratch / (stage + '.strict.ll'), adapter_headers, sealed[stage]['scope'])
                log = (owned / name / 'stderr.log').read_text()
                match = re.search(r'^SIMD16 value ABI: ' + stage + r' calls=(\d+) tags=(\d+) masked=(\d+)$', log, re.M)
                if not match or tuple(map(int, match.groups()[:2])) != (
                    sealed[stage]['scope']['checked_pointer_calls'], sealed[stage]['scope']['checked_tag_banks']
                ):
                    raise RuntimeError(f'{stage} bank proof disagrees with the sealed scope')
                detail['normalized_private_masked_accesses'] = int(match[3])
                results[name]['value'] = detail
            results[name]['output'] = product(output)
    finally:
        if results:
            receipt.parent.mkdir(parents=True, exist_ok=True)
            write_json(receipt, {'scope': 'SIMD16 compiler value ABI and adapter signature check', 'source': revision,
                                 'driver': product(Path(__file__).resolve()), 'pass': product(VALUE_SOURCE),
                                 'strict_scope': product(SCOPE), 'adapter_source': product(adapter),
                                 'sdk': descriptor['archive'], 'results': results})


if __name__ == '__main__':
    main()
