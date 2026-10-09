#!/usr/bin/env python3
"""Offline six-stage candidate; owns only this ignored out directory."""
import difflib
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import struct
import subprocess
import sys
import time
sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ORIGINAL = HERE/'original'
COPIED = HERE/'stage-src/sirius/kernels'
MODULES = HERE/'stage-modules'
PROGRAMS = HERE/'stage-programs'
PRESET = ROOT/'bin/linux-gcc/src/sirius/backend/retained'
for directory in [ORIGINAL, COPIED, MODULES, PROGRAMS]:
    directory.mkdir(parents=True, exist_ok=True)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run(args):
    result = subprocess.run(list(map(str, args)), cwd=HERE, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f'failed command {args}:\n{result.stdout}\n{result.stderr}')
    return result.stdout

def save(report):
    (HERE/'stage-artifact-report.json').write_text(json.dumps(report, indent=2) + '\n')

stages = [('Camera', 'camera', 'build_camera_program'),
          ('Transport', 'transport', 'build_transport_program'),
          ('Endpoint', 'endpoint', 'build_endpoint_program'),
          ('Dense', 'dense', 'build_dense_program'),
          ('Initialize', 'initialize', 'build_initialize_program'),
          ('RayCamera', 'ray_camera', 'build_ray_camera_program')]
original_paths = list(sorted((ROOT/'src/sirius/kernels').glob('retained*.slang')))
original_paths += [ROOT/'src/sirius/kernels/portable_binary32.h',
                   ROOT/'src/sirius/kernels/retained_program.py',
                   ROOT/'scripts/build-retained-kernels.py',
                   ROOT/'src/sirius/backend/CMakeLists.txt']
source_hashes = {str(p.relative_to(ROOT)): sha(p) for p in original_paths}
for path in original_paths:
    relative = path.relative_to(ROOT)
    snapshot = ORIGINAL/relative
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    snapshot.write_bytes(path.read_bytes())
    assert sha(snapshot) == source_hashes[str(relative)]
    if path.parent == ROOT/'src/sirius/kernels':
        (COPIED/path.name).write_bytes(snapshot.read_bytes())
helper = ROOT/'out/retained-qualification/portable-fp64-scalar/fp64_scalar.h'
assert sha(helper) == '5a348050a1b87bcd5b189de51240eac5753ba5d854919e7d35ce2d2499b7ff47'
(COPIED/'fp64_scalar.h').write_bytes(helper.read_bytes())
product_helper = HERE/'fp64_product.h'
assert sha(product_helper) == 'fda6079a13a45e2197a14ebae7286aa40da6f4dc1c8d6046a4e09d64e9d494b1'
(COPIED/'fp64_product.h').write_bytes(product_helper.read_bytes())
scalar = COPIED/'retained_scalar.slang'
old_scalar = scalar.read_text()
new_scalar = old_scalar.replace('#include "portable_binary32.h"', '#include "portable_binary32.h"\n#include "fp64_product.h"', 1)
for function, before, after in [('OrderedAdd32', 'PB32Add', 'PB64Add'),
                                ('OrderedSubtract32', 'PB32Subtract', 'PB64Subtract'),
                                ('OrderedMultiply32', 'PB32Multiply', 'PB64Multiply')]:
    old = f'RPScalar {function}(RPScalar a, RPScalar b) {{ return {before}(a, b); }}'
    new = f'RPScalar {function}(RPScalar a, RPScalar b) {{ return {after}(a, b); }}'
    assert new_scalar.count(old) == 1
    new_scalar = new_scalar.replace(old, new, 1)
scalar.write_text(new_scalar)
(HERE/'scalar-only.patch').write_text(''.join(difflib.unified_diff(
    old_scalar.splitlines(keepends=True), new_scalar.splitlines(keepends=True),
    fromfile='original/src/sirius/kernels/retained_scalar.slang', tofile='src/sirius/kernels/retained_scalar.slang')))
pair = COPIED/'retained_pair.slang'
old_pair = pair.read_text()
old_product = 'PB32Product product=PB32ProductResidual(RPBits(a),RPBits(b));'
assert old_pair.count(old_product) == 1
new_pair = old_pair.replace(old_product, 'PB32Product product=PB64ProductResidual(RPBits(a),RPBits(b));', 1)
old_comment = ("// The exact bounded integer product/residual avoids emulating Dekker's\n"
               " // repeated scalar operations. Both binary32 projections satisfy the same\n"
               " // five-minimum-subnormal enclosure used by the native product transforms.")
assert new_pair.count(old_comment) == 1
new_pair = new_pair.replace(old_comment,
    "// Isolated bounded binary64 product/residual under the reviewed provider\n"
    " // precision premise. Both projections preserve the unchanged\n"
    " // five-minimum-subnormal enclosure contract.", 1)
pair.write_text(new_pair)
(HERE/'product-only.patch').write_text(''.join(difflib.unified_diff(old_pair.splitlines(keepends=True), new_pair.splitlines(keepends=True), fromfile='original/src/sirius/kernels/retained_pair.slang', tofile='stage-src/sirius/kernels/retained_pair.slang')))
assert sha(COPIED/'portable_binary32.h') == '2f835dc2cbda7a62425977b31dc89155bc0393284d5ce3ace90f42c71bd076f2'
unchanged = {}
for path in original_paths:
    if path.parent == ROOT/'src/sirius/kernels' and path.name not in ('retained_scalar.slang', 'retained_pair.slang'):
        assert sha(COPIED/path.name) == source_hashes[str(path.relative_to(ROOT))]
        unchanged[path.name] = sha(COPIED/path.name)

spec = importlib.util.spec_from_file_location('candidate_retained_program', COPIED/'retained_program.py')
program_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(program_module)
# Preset artifacts are read-only provenance controls, never compiler outputs.
baseline_paths = [PRESET/f'retained_{stem}{suffix}.spv' for _, stem, _ in stages
                  for suffix in ['_portable', '_portable_fp64']]
baseline_paths += [PRESET/'retained_kernels.h']
baseline_hashes = {str(path.relative_to(ROOT)): sha(path) for path in baseline_paths}
baseline_header = (PRESET/'retained_kernels.h').read_text()
report = {
    'scope': 'offline finite feasibility candidate; no GPU, production integration or scientific qualification',
    'status': 'building', 'source_head': run(['git', '-C', ROOT, 'rev-parse', 'HEAD']).strip(),
    'source_branch': run(['git', '-C', ROOT, 'branch', '--show-current']).strip(),
    'original_sources_sha256': source_hashes, 'unchanged_copied_kernel_sources_sha256': unchanged,
    'reviewed_helper': {'original_path': str(helper.relative_to(ROOT)), 'sha256': sha(helper)},
    'modified_scalar_sha256': sha(scalar), 'modified_pair_sha256': sha(pair), 'product_patch_sha256': sha(HERE/'product-only.patch'), 'reviewed_product_helper_sha256': sha(product_helper), 'scalar_patch_sha256': sha(HERE/'scalar-only.patch'),
    'script_sha256': sha(Path(__file__)), 'preset_read_only_sha256': baseline_hashes,
    'compiler': {'path': '/opt/slang/bin/slangc', 'sha256': sha(Path('/opt/slang/bin/slangc')),
                 'version': run(['/opt/slang/bin/slangc', '-version']).strip(),
                 'optimization': '-O0 (unchanged canonical retained compiler policy)'},
    'artifact_policy': 'canonical compiler/profile/workgroup/program definitions; explicit binary64 RTE64 and NoContraction validation replace the pure-integer capability guard for this isolated candidate only',
    'device_calls': 0, 'dispatches': 0, 'stages': [],
}
save(report)

for kind, stem, builder in stages:
    program = getattr(program_module, builder)(parallel=True)
    terms = 4 if kind in ('Camera', 'RayCamera') else 5
    layers = len(program['layer_offsets']) - 1
    prefix = program.get('prefix_instructions', 0)
    shared_bytes = program['registers'] * terms * 4 + 4
    assert shared_bytes + 4 <= 16384
    program_words = [program['instructions'], program['registers']]
    if kind not in ('Camera', 'RayCamera'):
        program_words.append(len(program['outputs']))
    program_words += program['outputs'] + program['operations'] + program['layer_offsets']
    header_array = re.search(r'k' + kind + r'Program\{\{(.*?)\}\};', baseline_header, re.S)[1]
    preset_words = list(map(int, re.findall(r'(\d+)u', header_array)))
    assert program_words == preset_words, f'{kind} generated program differs from the canonical preset'
    program_path = PROGRAMS/f'retained_{stem}.bin'
    program_path.write_bytes(struct.pack('<' + str(len(program_words)) + 'I', *program_words))
    (PROGRAMS/f'retained_{stem}.json').write_text(json.dumps(program, indent=2) + '\n')
    words_per_row = ((512 if kind == 'Camera' else 576) + 4 * program['registers']
                     if kind in ('Camera', 'RayCamera') else
                     {'Transport': 2404, 'Endpoint': 769, 'Dense': 204, 'Initialize': 204}[kind] + 5 * program['registers'])
    preset_row_words = int(re.search(r'k' + kind + r'RowWords = (\d+);', baseline_header)[1])
    assert words_per_row == preset_row_words
    record = {'kind': kind, 'instructions': program['instructions'], 'registers': program['registers'],
              'layers': layers, 'terms': terms, 'prefix_instructions': prefix, 'row_words': words_per_row,
              'shared_memory_bytes': shared_bytes, 'program_sha256': sha(program_path),
              'program_matches_preset': True, 'modes': []}
    for wide in [False, True]:
        label = 'portable_fp64_products' if wide else 'portable_fp32_products'
        filename = f'retained_{stem}_{label}'
        raw = MODULES/(filename + '.compiler.spv')
        raw_assembly = MODULES/(filename + '.compiler.spvasm')
        patched_assembly = MODULES/(filename + '.patched.spvasm')
        output = MODULES/(filename + '.spv')
        disassembly = MODULES/(filename + '.spvasm')
        definitions = ['-DSIRIUS_RETAINED_PORTABLE=1']
        if wide:
            definitions.append('-DSIRIUS_RETAINED_FP64=1')
        definitions += [f'-DSIRIUS_RETAINED_REGISTERS={program["registers"]}',
                        f'-DSIRIUS_RETAINED_TERMS={terms}', '-DSIRIUS_RETAINED_LANES=64',
                        f'-DSIRIUS_RETAINED_LAYERS={layers}', f'-DSIRIUS_RETAINED_PREFIX={prefix}']
        command = ['/opt/slang/bin/slangc', COPIED/f'retained_{stem}.slang', *definitions,
                   '-I', COPIED, '-I', HERE/'stage-src', '-O0', '-target', 'spirv', '-profile', 'spirv_1_5',
                   '-entry', 'ComputeMain', '-stage', 'compute', '-o', raw]
        started = time.monotonic()
        run(command)
        run(['spirv-dis', raw, '-o', raw_assembly])
        text = raw_assembly.read_text()
        assert 'OpCapability RoundingModeRTE' not in text
        entry = re.search(r'OpEntryPoint GLCompute (%\S+)', text)[1]
        local = re.search(r'^.*OpExecutionMode ' + re.escape(entry) + r' LocalSize.*$', text, re.M)[0]
        assert local.strip() == f'OpExecutionMode {entry} LocalSize 64 1 1'
        text = text.replace('OpCapability Shader', 'OpCapability Shader\n               OpCapability RoundingModeRTE', 1)
        text = text.replace(local, local + '\n               OpExecutionMode ' + entry + ' RoundingModeRTE 64', 1)
        patched_assembly.write_text(text)
        run(['spirv-as', '--target-env', 'spv1.5', patched_assembly, '-o', output])
        run(['spirv-val', '--target-env', 'vulkan1.2', output])
        run(['spirv-dis', output, '-o', disassembly])
        # Validate exact final output, independently of intermediate dumps.
        text = disassembly.read_text()
        capabilities = re.findall(r'OpCapability (\S+)', text)
        assert len(capabilities) == 3 and set(capabilities) == {'Shader', 'Float64', 'RoundingModeRTE'}
        assert set(re.findall(r'OpTypeInt (\d+) [01]', text)) == {'32'}
        assert re.findall(r'OpTypeFloat (\d+)', text) == ['64']
        assert 'OpFConvert' not in text and ' Fma ' not in text
        double_type = re.search(r'(%\S+) = OpTypeFloat 64', text)[1]
        float_ops = re.findall(r'(%\S+) = (OpF(?!unction)\w+) (%\S+)', text)
        assert float_ops and {op for _, op, _ in float_ops} == {'OpFAdd', 'OpFSub', 'OpFMul'}
        assert all(t == double_type for _, _, t in float_ops)
        decorated = set(re.findall(r'OpDecorate (%\S+) NoContraction', text))
        assert all(result in decorated for result, _, _ in float_ops)
        entry = re.search(r'OpEntryPoint GLCompute (%\S+)', text)[1]
        assert set(re.findall(r'OpExecutionMode\S* (%\S+) (.*)', text)) == {
            (entry, 'LocalSize 64 1 1'), (entry, 'RoundingModeRTE 64')}
        assert set(re.findall(r'OpDecorate %\S+ DescriptorSet (\d+)', text)) == {'0'}
        assert sorted(re.findall(r'OpDecorate %\S+ Binding (\d+)', text)) == ['0', '1']
        assert set(re.findall(r'OpDecorate %\S+ ArrayStride (\d+)', text)) == {'4'}
        # Exact scratch words and uint32 status retain the canonical 16KiB bound.
        constants = {name: int(value) for name, value in re.findall(r'(%\S+) = OpConstant %(?:uint|int) (\d+)', text)}
        scratch_types = re.findall(r'(%\S+) = OpTypeArray %uint (%\S+)', text)
        assert any(constants.get(count) == program['registers'] * terms for _, count in scratch_types)
        mode = {'label': label, 'wide_product_definition': wide, 'path': str(output.relative_to(HERE)),
                'sha256': sha(output), 'bytes': output.stat().st_size,
                'disassembly_sha256': sha(disassembly), 'command': list(map(str, command)),
                'capabilities': capabilities, 'integer_widths': [32], 'float_widths': [64],
                'execution_modes': ['LocalSize64x1x1', 'RoundingModeRTE64'],
                'native_double_operation_count': len(float_ops), 'all_NoContraction': True,
                'no_Float32_or_Int64': True, 'validation': 'spirv-val --target-env vulkan1.2 PASS',
                'preparation_wall_seconds': time.monotonic() - started}
        record['modes'].append(mode)
        for temporary in [raw, raw_assembly, patched_assembly]:
            temporary.unlink()
        print(json.dumps({'stage': kind, 'mode': label, 'bytes': mode['bytes'], 'sha256': mode['sha256'],
                          'double_operations': mode['native_double_operation_count'], 'validated': True}), flush=True)
    record['public_mode_modules_identical'] = record['modes'][0]['sha256'] == record['modes'][1]['sha256']
    report['stages'].append(record)
    save(report)

report['original_sources_unchanged'] = all(sha(ROOT/name) == value for name, value in source_hashes.items())
report['preset_artifacts_unchanged'] = all(sha(ROOT/name) == value for name, value in baseline_hashes.items())
report['copied_sources_sha256'] = {str(p.relative_to(HERE)): sha(p) for p in sorted(COPIED.iterdir()) if p.is_file()}
assert report['original_sources_unchanged'] and report['preset_artifacts_unchanged']
report['status'] = 'complete'
report['pass'] = len(report['stages']) == 6 and all(len(v['modes']) == 2 for v in report['stages'])
save(report)
print(json.dumps({'complete': report['pass'], 'stages': len(report['stages']), 'modules': 12,
                  'device_calls': 0, 'dispatches': 0,
                  'all_public_modes_identical': all(v['public_mode_modules_identical'] for v in report['stages'])}), flush=True)
