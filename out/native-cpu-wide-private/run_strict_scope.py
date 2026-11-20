"""Verify the strict scope for the three already regenerated SIMD16 modules."""

from collections import Counter
from pathlib import Path, PurePosixPath
import json
import re
import resource
import shlex
import shutil
import subprocess
import sys
import tarfile

from resume_raw import BASELINE, GUARD, R, check_space, ir_inventory, pin, read_json, write_json


HISTORICAL_SDK = R / 'attestations/qualification/97b3a761/native-cpu-strict-ir-linking-compiler-gate/MANIFEST.json'
SOURCE = Path(__file__).with_name('strict_scope.cpp').resolve()
STAGES = ('transport', 'endpoint', 'dopri_phase')
TOOLS = ('clang++', 'llvm-config', 'opt')


def extract_sdk(archive, destination, descriptor):
    expected = descriptor['archive']
    if pin(archive) != {key: expected[key] for key in ('bytes', 'sha256')}:
        raise RuntimeError('The exact LLVM 23.1.2 SDK archive is required')
    members = {Path(descriptor['tools'][tool]['resolved']).relative_to(
        'bin/linux-gcc/native-cpu-strict-ir-linking/toolchain') for tool in TOOLS}
    members.add(Path('lib/libLLVMSupport.a'))
    # Read the xz via the native decoder and write only LLVM headers, the
    # three executables and the support archive containing the SHA member.
    includes = (Path('include/llvm'), Path('include/llvm-c'), Path('lib/clang/23/include'))
    with subprocess.Popen(['/usr/bin/xz', '-T4', '-dc', str(archive)], stdout=subprocess.PIPE) as decoder:
        with tarfile.open(mode='r|', fileobj=decoder.stdout) as bundle:
            for member in bundle:
                path = PurePosixPath(member.name)
                if path.is_absolute() or '..' in path.parts or path.parts[0] != descriptor['extracted_root']:
                    raise RuntimeError(f'Unexpected SDK archive path {path}')
                relative = Path(*path.parts[1:])
                if not member.isfile() or (relative not in members and not any(relative.is_relative_to(p) for p in includes)):
                    continue
                target = destination / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                with bundle.extractfile(member) as source, target.open('xb') as output:
                    shutil.copyfileobj(source, output, length=1024 * 1024)
                target.chmod(member.mode & 0o777)
        # The tar reader can stop at the terminator before xz has drained.
        while decoder.stdout.read(1024 * 1024):
            pass
        if decoder.wait() != 0:
            raise RuntimeError('LLVM SDK archive decompression failed')
    tools = {}
    for name in TOOLS:
        record = descriptor['tools'][name]
        resolved = Path(record['resolved']).relative_to('bin/linux-gcc/native-cpu-strict-ir-linking/toolchain')
        path = destination / resolved
        if pin(path) != {key: record[key] for key in ('bytes', 'sha256')}:
            raise RuntimeError(f'LLVM SDK tool fingerprint changed: {name}')
        tools[name] = path
    return tools


def build_action(scratch, tools):
    support = scratch / 'sdk/lib/libLLVMSupport.a'
    member = scratch / 'SHA256.cpp.o'
    with member.open('xb') as stream:
        subprocess.run(['/usr/bin/ar', 'p', str(support), member.name], stdout=stream, check=True)
    sha_archive = scratch / 'sha256-only.a'
    subprocess.run(['/usr/bin/ar', 'cr', str(sha_archive), str(member)], check=True)
    flags = shlex.split(subprocess.check_output([str(tools['llvm-config']), '--cxxflags'], text=True))
    return [str(tools['clang++']), '--driver-mode=g++', *flags, '-O2', '-fPIC', '-shared', str(SOURCE),
            str(sha_archive), '-Wl,--exclude-libs,ALL', '-o', str(scratch / 'sirius-wide-strict.so')]


def audit(raw, path):
    text = path.read_text()
    headers = re.findall(r'^define .*$', text, re.M)
    attrs = dict(re.findall(r'^attributes #(\d+) = (.*)$', text, re.M))
    if len(headers) != raw['definitions'] or any(
        not any('strictfp' in attrs[index].split() for index in re.findall(r'#(\d+)\b', header))
        for header in headers
    ):
        raise RuntimeError('Strictfp is missing on a SIMD16 definition')
    ops = Counter(re.findall(
        r'^\s+%[^=\n]+ = (fadd|fsub|fmul|fdiv|frem|fneg|fcmp|sitofp|uitofp|fptosi|fptoui|fpext|fptrunc)\b',
        text, re.M))
    if ops != Counter({'fneg': raw['fp_operations'].get('fneg', 0)}):
        raise RuntimeError(f'Unexpected ordinary FP operations {ops}')

    prefix = 'llvm.experimental.constrained.'
    expect = Counter({prefix + name + '.v16f32': raw['fp_operations'].get(name, 0)
                      for name in ('fadd', 'fsub', 'fmul')})
    predicates = raw['comparison_predicates']
    quiet = sum(predicates.get(name, 0) for name in ('oeq', 'une', 'uno'))
    expect[prefix + 'fcmp.v16f32'] = quiet
    expect[prefix + 'fcmps.v16f32'] = sum(predicates.values()) - quiet
    expect[prefix + 'sitofp.v16f32.v16i32'] = raw['fp_operations'].get('sitofp', 0)
    calls = Counter(re.findall(r'\bcall\s+[^\n]*?@([\w.$-]+)\(', text))
    constrained = Counter({name: count for name, count in calls.items() if name.startswith(prefix)})
    if constrained != +expect:
        raise RuntimeError(f'Constrained instruction inventory changed: {constrained}')
    classified = Counter()
    for line in re.findall(r'^.*\bcall\b[^\n]*?@llvm\.experimental\.constrained\.[^\n]+$', text, re.M):
        if 'metadata !"fpexcept.strict"' not in line or (
                '@llvm.experimental.constrained.fcmp' not in line and 'metadata !"round.dynamic"' not in line):
            raise RuntimeError('Constrained rounding or exception mode changed')
        cmp = re.search(r'@llvm\.experimental\.constrained\.(fcmp[s]?)\.v16f32\(.*metadata !"(\w+)",', line)
        if cmp:
            classified[cmp[1], cmp[2]] += 1
    if classified != Counter({('fcmp' if name in ('oeq', 'une', 'uno') else 'fcmps', name): count
                              for name, count in predicates.items()}):
        raise RuntimeError('Quiet/signaling classification changed')
    remaining = Counter({name: count for name, count in calls.items() if not name.startswith(prefix)})
    expected_calls = Counter({name.replace('sirius_wide_', 'sirius_wide_strict_', 1)
                              if name.startswith('sirius_wide_') else name: count
                              for name, count in raw['calls'].items()})
    actual_context, = re.findall(r'^%v16_varying_KernelContext_0 = type .+$', text, re.M)
    if remaining != expected_calls or actual_context != raw['context']:
        raise RuntimeError('Raw non-FP call inventory or context layout changed')
    return {'output': {'path': str(path.relative_to(R)), **pin(path)}, 'definitions': len(headers),
            'constrained_operations': dict(constrained), 'retained_fneg': ops.get('fneg', 0),
            'non_constrained_call_inventory_equal': True, 'raw_context_layout_equal': True}


def main():
    scratch, receipt = (Path(arg).resolve() for arg in sys.argv[1:])
    if not scratch.is_relative_to(R / 'out') or not receipt.is_relative_to(R / 'attestations') or receipt.exists():
        raise RuntimeError('Use a task workspace under out/ and a new manifest path under attestations/')
    if Path.cwd() != R or subprocess.check_output(['git', 'status', '--porcelain']).strip():
        raise RuntimeError('Run from the clean Sirius checkout')
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    check_space()
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    baseline = read_json(BASELINE)
    raw_inputs = {stage: R / 'bin/linux-gcc/native-cpu-wide-private' / (stage + '.ll') for stage in STAGES}
    inventories = {stage: ir_inventory(path) for stage, path in raw_inputs.items()}
    for stage, raw in inventories.items():
        expected = baseline['historical_transport'] if stage == 'transport' else baseline['results']['emit-' + stage]['product']
        if not all(raw[key] == value for key, value in expected.items()):
            raise RuntimeError(f'Ordinary IR pin differs for {stage}')
    descriptor = read_json(HISTORICAL_SDK)['setup']['official_sdk']
    tools = extract_sdk(scratch / 'LLVM-23.1.2-Linux-X64.tar.xz', scratch / 'sdk', descriptor)

    owned = scratch / 'scope'
    owned.mkdir()
    shutil.copyfile(GUARD / 'run_owned.py', owned / 'run_owned.py')
    (owned / 'owned_guard.py').symlink_to(GUARD / 'owned_guard.py')
    commands = {'build': build_action(scratch, tools)}
    for stage in STAGES:
        commands[stage] = [str(tools['opt']), '-load-pass-plugin=' + commands['build'][-1],
                           '-passes=sirius-wide-strict,verify', '-verify-each', '-S', str(raw_inputs[stage]),
                           '-o', str(scratch / (stage + '.strict.ll'))]
    actions = {}
    owner_inputs = [str(path) for path in (BASELINE, HISTORICAL_SDK, SOURCE, Path(__file__).resolve(), *tools.values())]
    headers = sorted(str(path) for directory in ('sdk/include', 'sdk/lib/clang/23/include')
                     for path in (scratch / directory).rglob('*') if path.is_file())
    for stage, command in commands.items():
        inputs = [commands['build'][-4], *headers] if stage == 'build' else [commands['build'][-1], str(raw_inputs[stage])]
        actions[stage] = dict(revision=revision, kind='diagnostic', expected_exit=0, software_provider=False,
                              limits=baseline['replay']['limits']['emit-' + (stage if stage != 'build' else 'transport')],
                              argv=command, frozen_inputs=owner_inputs + inputs,
                              scope='Pinned SIMD16 strict compiler scope')
    write_json(owned / 'actions.json', actions)
    results = {}
    try:
        for stage in commands:
            check_space()
            status = subprocess.run([sys.executable, '-B', str(owned / 'run_owned.py'), stage])
            owner_path = owned / stage / 'owner.json'
            if owner_path.exists():
                owner = read_json(owner_path)
                results[stage] = {key: owner[key] for key in (
                    'passed', 'source', 'wall_seconds', 'child_exit', 'stop_reason', 'whole_inputs_unchanged',
                    'peak_sampled_rss_kib') if key in owner}
                if not owner['passed']:
                    results[stage]['diagnostic'] = (owned / stage / 'stderr.log').read_text()[-2000:]
            status.check_returncode()
            output = Path(commands[stage][-1])
            if stage != 'build':
                raw = inventories[stage]
                results[stage]['raw'] = {key: raw[key] for key in ('path', 'bytes', 'sha256')}
                result = audit(raw, output)
                log = (owned / stage / 'stderr.log').read_text()
                observed = re.search(r'^SIMD16 strict scope: ' + stage +
                                     r' converted=(\d+) pointer_calls=(\d+) tagged_calls=(\d+) fmaxf=(\d+)$',
                                     log, re.M)
                pointer_calls = sum(count for name, count in raw['calls'].items() if name.startswith('sirius_wide_pointer_'))
                tagged_calls = sum(count for name, count in raw['calls'].items() if name.startswith('sirius_wide_pointer_tagged_'))
                if not observed or tuple(map(int, observed.groups())) != (
                    sum(result['constrained_operations'].values()), pointer_calls, tagged_calls,
                    raw['calls'].get('sirius_simd_fmaxf', 0),
                ):
                    raise RuntimeError('Pass report disagrees with the emitted IR')
                result.update(checked_pointer_calls=pointer_calls, checked_tag_banks=tagged_calls)
                results[stage]['scope'] = result
            else:
                results[stage]['plugin'] = pin(output)
    finally:
        if results:
            receipt.parent.mkdir(parents=True, exist_ok=True)
            write_json(receipt, {'scope': 'Strict scope compiler check only', 'source': revision,
                                 'driver': pin(Path(__file__)), 'pass': pin(SOURCE),
                                 'sdk': descriptor['archive'], 'results': results})


if __name__ == '__main__':
    main()
