"""Derive only same-eight-lane optional raw packets and a baseline packet selector."""
from pathlib import Path
import hashlib
import json
import re
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-v4-raw-packets'
PINS = {'out/native-cpu-local-identities/ispc-dense-compile-result.json': 'c1418cad5d85447280582b0fd4f81c735a075c210ebff8e0984ad3fd11accd56', 'out/native-cpu-local-identities/ispc-derivation.json': '0ab87b4d581c40981dc594b42528e39adcae25b5e7dc8dbe1c9483f13bf8d668', 'out/native-cpu-local-identities/derive.py': '095127a4b017bb6f3b86f603b49d4bdc06cbb5e5932f16fd82fa26856b85110f', 'out/native-cpu-local-identities/identity_store.isph': '0d464b6ce820abc38301d8f305ba6802ae04e5c13eff80c953bad118b786be2c', 'out/native-cpu-value-fact-stores/value_facts.isph': '7a245fb084f620a394f7920723b55d66b6c5986a21af78c78ea92cfb149c9f23', 'out/native-cpu-simd-packet-policy/rows_adapter.cpp': '57747ee37daa040d25b1e6f9411affb07fb399ba0003c921aca313810a10a581', 'out/native-cpu-simd-packet-policy/rows_bridge.h': 'd52b7cb2080ee35f0f06ea0df9e5e71f8694a3cc45024b486b5cb7a50b9dc282', 'out/native-cpu-simd-packet-policy/adapter-compile-result.json': 'd2c0d6bce1b3dc7ad39c1f004dda5f350efef64a062a041f9e87c969eaab5a0a', 'out/native-cpu-value-fact-stores/prepared-cpp.json': '0631a88a4f7deb207073ee2c4219526ef77f2b668e128d6587cae24226f5c053', 'out/native-cpu-endpoint-header-provider/provider/source-contract.json': '468f9ce2b6b7d2c77327f3a8f3bd9c996855f78e489fb83b0b0c6bd640db57bf', 'out/native-cpu-v4-provider/provider/source-contract.json': '02b60b2e49a084a326d31a5843611b84607b0c4fb650f00cb968d7b4cf68a121', 'out/native-cpu-simd-rows/compiler-environment.json': 'f225ac34c573c886ee0163f25c511ed97029de5e6f0629f49d1001524bd861eb', 'out/native-cpu-simd-rows/compiler-acquisition.json': '1af00eea3a09ff7772324105fdbdcc001d62491a22c2a68ce741f3710a49e20b', 'bin/linux-gcc/native-cpu-simd-rows/compiler/ispc.html': '11ca3d2a0584083270fb77cb2dd95ca6d048ec0f5b3ac25a9638935a2c9652fc', 'out/native-cpu-v4-arithmetic/optional.cpp': 'bcd7804fc7935de78846bf6317d9545b7df78afa934957d74691e74927b028d5', 'out/native-cpu-v4-arithmetic/dispatch.cpp': 'ab1512b576709665846a60b5e68f0c4dc9dbe86c08a84e794c18622037117445'}
STAGES = ('transport', 'endpoint', 'dopri_phase')
FLAGS = ['-O2', '--target=avx2-i32x8', '--cpu=core-avx2', '--arch=x86-64',
         '--addressing=32', '--PIC', '--opt=disable-fma', '--math-lib=system']


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def key(path):
    return str(path.relative_to(R)) if path.is_relative_to(R) else str(path)


def checked(row):
    path = R / row['path']
    assert pin(path) == {name: row[name] for name in ('bytes', 'sha256')}, path


def derive():
    assert W == R / 'out/native-cpu-v4-raw-packets'
    assert not (W / 'derivation.json').exists() and not (W / 'rows_adapter.cpp').exists()
    assert not B.exists()
    for name, expected in PINS.items():
        assert pin(R / name)['sha256'] == expected, name
    raw = doc(R / 'out/native-cpu-local-identities/ispc-dense-compile-result.json')
    original = doc(R / 'out/native-cpu-local-identities/ispc-derivation.json')
    ap = doc(R / 'out/native-cpu-endpoint-header-provider/provider/source-contract.json')
    ao = doc(R / 'out/native-cpu-v4-provider/provider/source-contract.json')
    adapter = doc(R / 'out/native-cpu-simd-packet-policy/adapter-compile-result.json')
    assert raw['all_objects_compiled'] and raw['all_prepared_source_bytes_unchanged']
    assert raw['new_raw_compilations'] == 3 and original['instructions'] == 9380
    assert ap['native_objects'][4:7] == [row['object'] for row in raw['objects']]
    assert ap['adapter_input_objects'] == ao['adapter_input_objects']
    inputs = dict(raw['source_inputs'])
    for name, expected in ap['source_inputs'].items():
        assert name not in inputs or inputs[name] == expected, name
        inputs[name] = expected
    for name in PINS:
        inputs[name] = pin(R / name)
    inputs[key(Path(__file__).resolve())] = pin(Path(__file__).resolve())
    parts = ap['adapter_input_objects']
    math = [row for row in parts if Path(row['path']).name in ('optional.o', 'v3.o', 'dispatch.o')]
    assert len(math) == 3
    reused = ap['native_objects'][:7] + math
    for row in reused + [parts[0]]:
        checked(row)
        inputs[row['path']] = {name: row[name] for name in ('bytes', 'sha256')}
    compiler = Path(raw['objects'][0]['argv'][0])
    assert pin(compiler)['sha256'] == '1813a6b83442dfcc23a1a82a641f22452fedbff7b429877e2f31957957e0168e'
    guide = R / 'bin/linux-gcc/native-cpu-simd-rows/compiler/ispc.html'
    text = guide.read_text()
    assert 'avx512skx-x8' in text and 'disable-fma' in text and 'disable-zmm' in text
    names = ['sirius_native_' + kind + '_' + op + '8' for kind in ('dense', 'tagged')
             for op in ('add', 'multiply', 'subtract', 'divide', 'sqrt')]
    mapping = {name: name.replace('sirius_native_', 'sirius_native_v4_', 1) for name in names}
    sources, derived, commands, stage_proof = {}, {}, [], {}
    literal_total = 0
    for slug in STAGES:
        receipt, = [row for row in raw['objects'] if row['stage'] == slug]
        before_path = R / original['stages'][slug]['derived_path']
        before = before_path.read_bytes()
        assert pin(before_path) == original['stages'][slug]['derived']
        assert receipt['argv'][1] == str(before_path) and receipt['argv'][2:-3] == FLAGS
        assert receipt['argv'][-3:-1] == ['--emit-obj', '-o']
        old_export, new_export = 'sirius_simd_' + slug + '_packet', 'sirius_simd_v4_' + slug + '_packet'
        assert len(re.findall(rb'export void ' + old_export.encode() + rb'\(', before)) == 1
        for name in names:
            assert re.search(rb'extern "C" void ' + name.encode() + rb'\(', before), name
        symbols = {old_export: new_export, **mapping}
        prefix = ('// Only compile-time symbol bindings differ; original source below is exact.\n' +
                  ''.join('#define ' + a + ' ' + b + '\n' for a, b in symbols.items())).encode()
        after = prefix + before
        assert after[len(prefix):] == before
        calls = [line for line in before.splitlines() if line.strip().startswith(b'if(!EvaluateOneStatic_0<')]
        assert len(calls) == original['stages'][slug]['literal_instruction_calls_preserved']
        assert [line for line in after.splitlines() if line.strip().startswith(b'if(!EvaluateOneStatic_0<')] == calls
        literal_total += len(calls)
        destination = B / 'ispc-source' / (slug + '.ispc')
        sources[destination] = after
        argv = list(receipt['argv'])
        argv[1], argv[-1] = str(destination), str(B / 'objects' / (slug + '.o'))
        argv[argv.index('--target=avx2-i32x8')] = '--target=avx512skx-x8'
        argv[argv.index('--cpu=core-avx2')] = '--cpu=skx'
        assert argv[2:-3] == [flag.replace('avx2-i32x8', 'avx512skx-x8').replace('core-avx2', 'skx') for flag in FLAGS]
        commands.append(dict(stage=slug, argv=argv, output=argv[-1]))
        stage_proof[slug] = dict(basis_path=key(before_path), basis=pin(before_path),
            derived_path=key(destination), prefix=prefix.decode(), preprocessor_symbol_map=symbols,
            inverse_restores_exact_AN_source=True, literal_instruction_calls=len(calls),
            original_compile_argv=receipt['argv'],
            both_C_export_and_masked_ISPC_export_renamed_by_same_source_identifier=True)
        inputs[key(before_path)] = pin(before_path)
    assert literal_total == 9380
    basis = R / 'out/native-cpu-simd-packet-policy/rows_adapter.cpp'
    before = basis.read_text()
    after = before.replace('#include "rows_bridge.h"', '#include "../native-cpu-simd-packet-policy/rows_bridge.h"', 1)
    assert after != before
    marker = 'static_assert(sizeof(uint32_t) == 4);'
    declarations = 'extern "C" bool sirius_native_v4_supported() noexcept;\n'
    for slug in STAGES:
        declarations += ('extern "C" void sirius_simd_v4_' + slug + '_packet(uint32_t *input, uint32_t input_words,\n'
                         '    uint32_t *output, uint32_t output_words, uint32_t first_x, uint32_t groups);\n')
    assert after.count(marker) == 1
    after = after.replace(marker, declarations + '\n' + marker, 1)
    replacements = []
    for slug in STAGES:
        call = ('    sirius_simd_' + slug + '_packet(const_cast<uint32_t*>(static_cast<const uint32_t*>(input)),\n'
                '        input_bytes / 4U, static_cast<uint32_t*>(output), output_bytes / 4U, first_x, groups);')
        changed = ('    if (sirius_native_v4_supported())\n' + call.replace('    sirius_simd_', '        sirius_simd_v4_', 1) +
                   '\n    else\n' + '    ' + call)
        assert after.count(call) == 1
        after = after.replace(call, changed, 1)
        replacements.append((call, changed))
    inverse = after.replace(declarations + '\n', '', 1)
    for call, changed in replacements:
        assert inverse.count(changed) == 1
        inverse = inverse.replace(changed, call, 1)
    inverse = inverse.replace('#include "../native-cpu-simd-packet-policy/rows_bridge.h"', '#include "rows_bridge.h"', 1)
    assert inverse == before
    sources[W / 'rows_adapter.cpp'] = after.encode()
    adapter_row, = adapter['objects']
    assert adapter_row['exit'] == 0 and adapter['all_prepared_source_bytes_unchanged']
    assert adapter_row['object'] == parts[0]
    adapter_argv = list(adapter_row['argv'])
    assert '-march=x86-64' in adapter_argv and '-ffp-contract=off' in adapter_argv
    adapter_argv[adapter_argv.index('-c') + 1] = str(W / 'rows_adapter.cpp')
    adapter_argv[adapter_argv.index('-o') + 1] = str(B / 'objects/rows_adapter.o')
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    for path, data in sources.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            stream.write(data)
        derived[key(path)] = pin(path)
    environment = doc(R / 'out/native-cpu-simd-rows/compiler-environment.json')
    result = dict(prepared_only=True, compiled=False, executed=False, source_inputs=inputs,
        derived=derived, stages=stage_proof, literal_instruction_calls_preserved=literal_total,
        raw_compilations=commands, adapter_compile_argv=adapter_argv, reused_objects=reused,
        original_adapter=parts[0], original_adapter_compile_argv=adapter_row['argv'],
        adapter_inverse_restores_original=True, component_manifest=ap['component_manifest'],
        compiler=dict(path=key(compiler), **pin(compiler), version_output=environment['version_output']),
        local_guide=dict(path=key(guide), **pin(guide)),
        read_only_help_observation=dict(commands=['ispc --version', 'ispc --help', 'ispc --support-matrix'],
            supported_target='avx512skx-x8', supported_cpu='skx', Linux_x86_64_supported=True,
            eight_program_instances=True, disable_fma_supported=True,
            disable_zmm_is_not_a_replacement_for_eight_lane_target=True,
            provenance='Observed from this exact pinned compiler before derivation; derivation does not invoke a compiler.'),
        required_postcompile_gate=dict(before_any_link_or_runtime=True,
            inspect_all_executable_sections=True, reject_unknown_or_non_v4_instructions=True,
            reject_FMA_x87_FP_control_writes_and_reassociation=True,
            guard='Unchanged baseline-compiled sirius_native_v4_supported uses the sealed OS-aware x86-64-v4 builtin; original groups8/v3/admission ordering precedes it.',
            target_features=['AVX512F', 'AVX512DQ', 'AVX512CD', 'AVX512BW', 'AVX512VL'],
            cpu_feature_scope='skx tuning can enable instructions outside v4. Help support is not an emitted-ISA proof; ADX/cache-flush/random/MPX/TSX/PKU/SGX/XSAVE/monitor and every unknown extension must be rejected unless separately within the sealed guard.'),
        proof='Remove only each recorded #define prefix to recover all original raw source bytes. Every literal tuple, arithmetic expression, table, context, mask, programCount-dependent eight-lane stride and C pointer/u32 parameter remains. Only optional raw extern symbol bindings bypass repeated math dispatch after one baseline packet guard; original raw3 continue using their original dispatch functions. Adapter inverse restores original bytes; all scalar admission/tail/ordinary fallback paths and FP scopes remain.',
        scope='Prospective three raw AVX512SKX-x8 compiles plus one baseline adapter compile only. Exact original AO optional/v3/dispatcher math is reused; exact-product-classes math is excluded. Original AP scalar/admission and AN raw3 remain available. No source regeneration of shaders, sixteen-lane ABI, FMA, fast math, FP-state edits, scheduling/count/cache/governor/buffer/host/scene changes, compilation, linking, workloads, hit-rate/benefit/frame or scientific/release acceptance. Actual byte/ISA/relocation/ABI/mask/FP/whole-result qualification remains required.')
    with (W / 'derivation.json').open('x') as stream:
        stream.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(source_only=True, raw_sources=3, prospective_compilations=4,
                         literal_instruction_calls=literal_total, reused_objects=len(reused), executed=False)))


if __name__ == '__main__':
    assert sys.argv[1:] == ['derive']
    derive()
