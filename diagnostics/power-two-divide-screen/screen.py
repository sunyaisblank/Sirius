"""Bounded issue-99 compiler screen; production files are never edited."""
from pathlib import Path
import hashlib
import importlib.util
import json
import re
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
PRODUCTS = ROOT / 'bin/linux-gcc/power-two-divide-screen'
PRIOR = ROOT / 'attestations/software-vulkan/4bcd987/guarded-product-high-jit'
OLD = PRIOR / 'files/out/guarded-product-high-runtime'
REV = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
VARIANTS = [('endpoint', False), ('endpoint', True), ('transport', True),
            ('dopri_phase', True), ('dense', False)]

def record(path):
    path = Path(path)
    with path.open('rb') as f:
        digest = hashlib.file_digest(f, 'sha256').hexdigest()
    return {'bytes': path.stat().st_size, 'sha256': digest}

def write(name, value):
    (WORK / name).write_text(json.dumps(value, indent=2) + '\n')

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value

def prepare():
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == REV
    assert not subprocess.check_output(['git', 'status', '--porcelain=v1'])
    assert not PRODUCTS.exists()
    PRODUCTS.mkdir()
    (PRODUCTS / 'modules').mkdir()
    for owner in ['stock', 'candidate']:
        target = WORK / owner / 'src/sirius/kernels'
        target.mkdir(parents=True)
        for p in (ROOT / 'src/sirius/kernels').iterdir():
            if p.is_file() and p.suffix in ['.h', '.slang', '.py']:
                shutil.copyfile(p, target / p.name)
    header = WORK / 'candidate/src/sirius/kernels/portable_binary32.h'
    text = header.read_text()
    needle = '    int exponent = aa.exponent - bb.exponent;\n    uint remainder = aa.significand;'
    assert text.count(needle) == 1
    addition = ('    int exponent = aa.exponent - bb.exponent;\n'
                '    // Exact quotient: the original radix loop yields (aa.significand << 3, 0).\n'
                '    if (bb.significand == 0x00800000u)\n'
                '        return PB32Round(sign, exponent, aa.significand << 3);\n'
                '    uint remainder = aa.significand;')
    header.write_text(text.replace(needle, addition, 1))
    write('candidate-change.json', {'source_revision': REV, 'issue': 99,
        'original_needle': needle, 'replacement': addition,
        'scope': 'Only candidate PB32Divide changes; no production edit or adoption.'})
    (WORK / 'candidate.patch').write_text(subprocess.run(
        ['diff', '-u', str(WORK / 'stock/src/sirius/kernels/portable_binary32.h'), str(header)],
        capture_output=True, text=True).stdout)
    for name in ['run_owned.py', 'owned_guard.py', 'capture_jit.py', 'analyze_capture.py']:
        shutil.copyfile(OLD / name, WORK / name)
    capture = (OLD / 'run_capture.py').read_text().replace('guarded-product-high-runtime', 'power-two-divide-screen')
    (WORK / 'run_capture.py').write_text(capture)
    probe = PRIOR / 'files/bin/linux-gcc/guarded-product-high-runtime/short_input_jit_probe'
    expected = json.loads((OLD / 'capture-stock-endpoint-normal/inputs-before.json').read_text())
    assert record(probe) == {k: expected['bin/linux-gcc/guarded-product-high-runtime/short_input_jit_probe'][k] for k in ['bytes', 'sha256']}
    shutil.copy2(probe, PRODUCTS / 'short_input_jit_probe')
    (PRODUCTS / 'short_input_jit_probe').chmod(0o755)
    providers = ['/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so',
                 '/usr/lib/x86_64-linux-gnu/libLLVM.so.20.1',
                 '/usr/share/vulkan/icd.d/lvp_icd.json']
    for name in providers:
        assert record(name) == {k: expected[name][k] for k in ['bytes', 'sha256']}, name
    source_files = [p for owner in ['stock', 'candidate'] for p in (WORK / owner).rglob('*') if p.is_file()]
    scripts = [p for p in WORK.glob('*.py')]
    common = [str(p.resolve()) for p in source_files + scripts]
    common += providers + [str(ROOT / 'attestations/software-vulkan/b5d2abf/optimized-ir-stage-code-binding/optimizer_observer.so')]
    cache = (ROOT / 'bin/linux-gcc/CMakeCache.txt').read_text()
    tool = lambda n: re.search(r'^' + n + r':FILEPATH=(.*)$', cache, re.M)[1]
    tools = [tool(n) for n in ['SIRIUS_SLANGC', 'SIRIUS_SPIRV_AS', 'SIRIUS_SPIRV_DIS', 'SIRIUS_SPIRV_VAL', 'SIRIUS_SPIRV_OPT']]
    common += tools + [str(PRODUCTS / 'short_input_jit_probe')]
    actions = {'compile-modules': {'revision': REV, 'kind': 'diagnostic', 'expected_exit': 0,
        'limits': {'seconds': 300, 'rss_kib': 4 * 1024**2, 'reserve_kib': 2 * 1024**2},
        'argv': ['/usr/bin/python3', '-B', str(WORK / 'screen.py'), 'compile'],
        'frozen_inputs': common,
        'scope': 'Complete source-frozen stock/candidate SPIR-V with unchanged production builder and validation; no arithmetic or speed verdict.'}}
    for stage, normal in VARIANTS:
        variant = stage + ('-normal' if normal else '-core')
        actions['capture-candidate-' + variant] = {'revision': REV, 'kind': 'diagnostic', 'expected_exit': 0,
            'software_provider': True,
            'limits': {'seconds': 300, 'rss_kib': 4 * 1024**2, 'reserve_kib': 2 * 1024**2},
            'argv': ['/usr/bin/python3', '-B', str(WORK / 'run_capture.py'), '{folder}',
                     str(PRODUCTS / 'modules' / ('candidate-' + variant + '.spv'))],
            'frozen_inputs': common + [str(PRODUCTS / 'modules' / ('candidate-' + variant + '.spv'))],
            'scope': 'Short-input complete-module ordinary LLVM/JIT code binding only; guard returns before retained arithmetic; no numerical/performance qualification.'}
    write('actions.json', actions)
    write('source-freeze.json', {str(p.relative_to(WORK)): record(p) for p in source_files + scripts})
    write('preparation.json', {'source_revision': REV, 'issue': 99, 'production_unchanged': True,
        'stock_provider_and_probe_reused_exact': True, 'prior_evidence': str(PRIOR),
        'tools': {p: record(p) for p in tools}, 'source_files_per_owner': len(source_files) // 2})
    print('Prepared exact stock and one-branch candidate, original observer/probe and unchanged provider inputs.', flush=True)

def compile_modules():
    for name, expected in json.loads((WORK / 'source-freeze.json').read_text()).items():
        assert record(WORK / name) == expected
    builder = module('issue99_builder', ROOT / 'scripts/build-retained-kernels.py')
    programs = module('issue99_programs', WORK / 'stock/src/sirius/kernels/retained_program.py')
    cache = (ROOT / 'bin/linux-gcc/CMakeCache.txt').read_text()
    tool = lambda n: re.search(r'^' + n + r':FILEPATH=(.*)$', cache, re.M)[1]
    compiler, assembler, disassembler, validator, optimizer = [tool(n) for n in
        ['SIRIUS_SLANGC', 'SIRIUS_SPIRV_AS', 'SIRIUS_SPIRV_DIS', 'SIRIUS_SPIRV_VAL', 'SIRIUS_SPIRV_OPT']]
    results = []
    for stage, normal in VARIANTS:
        program = getattr(programs, 'build_' + stage + '_program')(parallel=True)
        kwargs = {}
        if stage in ['endpoint', 'transport']:
            special = getattr(programs, 'build_schwarzschild_' + stage + '_program')(parallel=True)
            encode = lambda p: [p['instructions'], p['registers'], len(p['outputs'])] + p['outputs'] + p['operations'] + p['layer_offsets']
            dimensions = (len(encode(program)), len(encode(special)), len(special['layer_offsets']) - 1)
            if stage == 'endpoint':
                dimensions += (special['prefix_instructions'],)
            kwargs[stage + '_specialized'] = dimensions
        variant = stage + ('-normal' if normal else '-core')
        texts = {}
        for owner in ['stock', 'candidate']:
            source = WORK / owner / 'src/sirius/kernels' / ('retained_' + stage + '.slang')
            output = PRODUCTS / 'modules' / (owner + '-' + variant + '.spv')
            builder.compile_shader(source, output, compiler, assembler, disassembler, validator,
                program['registers'], 5, len(program['layer_offsets']) - 1,
                program.get('prefix_instructions', 0), portable=True, optimizer=optimizer,
                normal_sum32=normal, **kwargs)
            assembly = output.with_suffix('.final.spvasm')
            subprocess.run([disassembler, str(output), '-o', str(assembly)], check=True)
            texts[owner] = assembly.read_text()
        canonical = ROOT / 'bin/linux-gcc/src/sirius/backend/retained' / ('retained_' + stage + ('_portable_normal_sum.spv' if normal else '_portable.spv'))
        stock = PRODUCTS / 'modules' / ('stock-' + variant + '.spv')
        candidate = PRODUCTS / 'modules' / ('candidate-' + variant + '.spv')
        assert stock.read_bytes() == canonical.read_bytes()
        assert candidate.read_bytes() != stock.read_bytes()
        shape = lambda text: sorted(re.findall(r'OpExecutionMode %\S+ (.*)|OpDecorate %\S+ ((?:Binding|DescriptorSet|ArrayStride) \d+)|OpMemberDecorate %\S+ (\d+ Offset \d+)', text))
        assert shape(texts['stock']) == shape(texts['candidate'])
        for instruction in ['OpControlBarrier', 'OpMemoryBarrier']:
            assert texts['stock'].count(instruction) == texts['candidate'].count(instruction)
        results.append({'variant': variant, 'stock': record(stock), 'candidate': record(candidate),
                        'production_policy_passed': True, 'stock_equals_current_canonical': True,
                        'mode_binding_stride_member_offset_inventory_exact': True})
        print(variant, 'compiled and structurally checked', flush=True)
    write('compiler-result.json', {'passed': True, 'modules': results,
        'scope': 'Full modules pass unchanged production structural policy; decoration inventories are supporting ABI evidence, not complete semantic equivalence.'})

if __name__ == '__main__':
    {'prepare': prepare, 'compile': compile_modules}[sys.argv[1]]()
