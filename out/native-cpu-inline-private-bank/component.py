"""Derive three inline opcode stages and zero-only private destination slots."""
from pathlib import Path
import hashlib
import json
import re
import sys
import zipfile

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-inline-private-bank'
STAGES = ('transport', 'endpoint', 'dopri_phase')
PARENT = R / 'attestations/qualification/5d6ba43c/native-cpu-packed-math-boundary-finite-compatibility/MANIFEST.json'
MATH = R / 'attestations/qualification/ca334e6e/native-cpu-exact-product-zero-residual-finite-compatibility/MANIFEST.json'
PINS = {
    str(PARENT.relative_to(R)): '12489ed18cdc6d429538d9c9b5b09ea685dbdbe7140cea83f3073a781647ca5d',
    str(MATH.relative_to(R)): '44d228597b02d724ad7a5b0ebfe8d69073bf306977bacd3a93013d1462213a66',
    'out/native-cpu-packed-math-provider/provider/linked-native-proof.json': 'e11cd5d2752722faf4385d7fa231adfcb28e9df3e02d83eab979937b57e79f7d',
    'out/native-cpu-packed-math-provider/provider/source-contract.json': '45949340a43c89f78194064b6a4cc45fff844bedc58ef96e9f9e630862e01ecf',
    'bin/linux-gcc/native-cpu-simd-rows/compiler/ispc.html': '11ca3d2a0584083270fb77cb2dd95ca6d048ec0f5b3ac25a9638935a2c9652fc',
    'out/native-cpu-packed-math-boundary/packinghelper.isph': 'b22326460ff2fdec4b3e507228f5ecd792cea82e33cf1f599647d893691cd874',
    'out/native-cpu-local-identities/identity_store.isph': '0d464b6ce820abc38301d8f305ba6802ae04e5c13eff80c953bad118b786be2c',
    'out/native-cpu-value-fact-stores/value_facts.isph': '7a245fb084f620a394f7920723b55d66b6c5986a21af78c78ea92cfb149c9f23',
}


def doc(path):
    return json.loads(path.read_text())


def data_pin(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def key(path):
    return str(path.relative_to(R)) if path.is_relative_to(R) else str(path)


def bundle(manifest):
    assert manifest['finite_component_passed'] and not manifest['provider_adopted']
    path = R / manifest['archive']['path']
    assert pin(path) == {name: manifest['archive'][name] for name in ('bytes', 'sha256')}
    return path


def derive():
    assert W == R / 'out/native-cpu-inline-private-bank'
    assert not (W / 'derivation.json').exists() and not (W / 'packinghelper.isph').exists()
    assert not B.exists()
    for name, sha in PINS.items():
        assert pin(R / name)['sha256'] == sha, name
    manifest, math = doc(PARENT), doc(MATH)
    archive, math_archive = bundle(manifest), bundle(math)
    inputs = {name: pin(R / name) for name in PINS}
    inputs.update({key(path): pin(path) for path in (archive, math_archive, Path(__file__).resolve())})
    static = doc(R / 'out/native-cpu-packed-math-provider/provider/linked-native-proof.json')
    provider = doc(R / 'out/native-cpu-packed-math-provider/provider/source-contract.json')
    assert static['passed'] and static['checked_optional_math_splitter_retention_dependencies'] == 36
    references = [row['object'] for row in static['raw_v4_objects']]
    assert len(references) == 3
    assert manifest['math_component_manifest'] == dict(path=key(MATH), **pin(MATH))
    optional = manifest['new_optional_math']
    assert optional in provider['adapter_input_objects']
    assert {name: optional[name] for name in ('bytes', 'sha256')} == math['products'][optional['path']]
    guide = R / 'bin/linux-gcc/native-cpu-simd-rows/compiler/ispc.html'
    assert 'inline</tt> is not a hint, but forces inlining.' in guide.read_text()
    sources, stages, commands = {}, {}, []
    with zipfile.ZipFile(archive) as zipped:
        def member(name):
            data = zipped.read(name)
            assert data_pin(data) == manifest['members'][name], name
            return data
        basis = json.loads(member('out/native-cpu-packed-math-boundary/derivation.json'))
        old_helper = R / 'out/native-cpu-packed-math-boundary/packinghelper.isph'
        helper_before = member(key(old_helper)).decode()
        assert helper_before.encode() == old_helper.read_bytes()
        payload_before = ('            packed[component*programCount+programIndex] =\n'
                          '                bank[(destination*5U+component)*programCount+programIndex];')
        payload_after = '            packed[component*programCount+programIndex] = 0U;'
        tag_before = '        packedTags[programIndex] = tags[destination*programCount+programIndex];'
        tag_after = '        packedTags[programIndex] = 0U;'
        helper_after = helper_before
        for before, after in ((payload_before, payload_after), (tag_before, tag_after)):
            assert helper_after.count(before) == 1
            helper_after = helper_after.replace(before, after, 1)
        inverse = helper_after.replace(payload_after, payload_before, 1).replace(tag_after, tag_before, 1)
        assert inverse == helper_before
        assert helper_before[helper_before.index('static inline void SiriusPackedMathCopyBack_0'): ] == helper_after[helper_after.index('static inline void SiriusPackedMathCopyBack_0'): ]
        helper_path = W / 'packinghelper.isph'
        sources[helper_path] = helper_after.encode()
        old_include, new_include = '#include "' + str(old_helper) + '"', '#include "' + str(helper_path) + '"'
        total = 0
        for slug, reference in zip(STAGES, references):
            source_path = R / basis['stages'][slug]['derived_path']
            before = member(key(source_path)).decode()
            assert pin(source_path) == basis['derived'][key(source_path)] == data_pin(before.encode())
            assert source_path.read_text() == before
            header = 'static bool EvaluateOneStatic_0('
            inline = 'static inline bool EvaluateOneStatic_0('
            assert before.count(header) == 1 and inline not in before
            assert before.count(old_include) == 1 and new_include not in before
            after = before.replace(header, inline, 1).replace(old_include, new_include, 1)
            assert after.replace(inline, header, 1).replace(new_include, old_include, 1) == before
            literals = lambda text: [line for line in text.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
            assert literals(before) == literals(after)
            count = len(literals(before))
            assert count == basis['stages'][slug]['literal_instruction_calls']
            total += count
            # Keep every layer declaration/call and every original source byte
            # except the sole opcode qualifier and required include-path binding.
            layer_declarations = re.findall(r'^static void EvaluateStaticLayer_[^\n]+', before, re.M)
            assert layer_declarations == re.findall(r'^static void EvaluateStaticLayer_[^\n]+', after, re.M)
            assert len(layer_declarations) in (159, 195, 45)
            original_scopes = [row for row in static['raw_v4_objects'][len(stages)]['functions']
                               if row['symbol'].startswith('EvaluateOneStatic_0')]
            assert len(original_scopes) == {'transport': 6, 'endpoint': 6, 'dopri_phase': 4}[slug]
            assert reference['path'] == 'bin/linux-gcc/native-cpu-packed-math-boundary/objects/' + slug + '.o'
            assert pin(R / reference['path']) == {name: reference[name] for name in ('bytes', 'sha256')}
            inputs[key(source_path)] = pin(source_path)
            inputs[reference['path']] = pin(R / reference['path'])
            destination = B / 'ispc-source' / (slug + '.ispc')
            sources[destination] = after.encode()
            original, = [row for row in basis['raw_compilations'] if row['stage'] == slug]
            argv = list(original['argv'])
            assert argv[1] == str(source_path)
            assert argv[2:-3] == ['-O2', '--target=avx512skx-x8', '--cpu=skx', '--arch=x86-64',
                                  '--addressing=32', '--PIC', '--opt=disable-fma', '--math-lib=system']
            assert argv[-3:-1] == ['--emit-obj', '-o']
            argv[1], argv[-1] = str(destination), str(B / 'objects' / (slug + '.o'))
            commands.append(dict(stage=slug, argv=argv, output=argv[-1]))
            stages[slug] = dict(basis_path=key(source_path), basis=pin(source_path),
                derived_path=key(destination), literal_instruction_calls=count,
                layer_scopes_unchanged=len(layer_declarations), whole_source_inverse_restores_parent=True,
                sole_opcode_inline_qualifier=True, helper_include_rebound=True,
                expected_removed_EvaluateOneStatic_functions=original_scopes,
                original_compile_argv=original['argv'],
                preprocessor_symbol_map=basis['stages'][slug]['preprocessor_symbol_map'])
        assert total == 9380
        compiler = R / basis['compiler']['path']
        assert pin(compiler) == {name: basis['compiler'][name] for name in ('bytes', 'sha256')}
        inputs[key(compiler)] = pin(compiler)
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    derived = {}
    for path, data in sources.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            stream.write(data)
        derived[key(path)] = pin(path)
    reference_paths = {row['path'] for row in references}
    reused = provider['native_objects'][:7] + [row for row in provider['adapter_input_objects']
                                               if row['path'] not in reference_paths]
    receipt = dict(prepared_only=True, compiled=False, executed=False, source_inputs=inputs,
        derived=derived, stages=stages, raw_compilations=commands, reference_objects=references,
        reused_objects=reused, literal_instruction_calls_preserved=total,
        compiler=basis['compiler'],
        parent_component_manifest=dict(path=key(PARENT), **pin(PARENT)),
        math_component_manifest=dict(path=key(MATH), **pin(MATH)),
        parent_static_proof=dict(path='out/native-cpu-packed-math-provider/provider/linked-native-proof.json',
            **pin(R / 'out/native-cpu-packed-math-provider/provider/linked-native-proof.json')),
        packing_helper=dict(basis_path=key(old_helper), basis=data_pin(helper_before.encode()),
            derived_path=key(helper_path), derived=data_pin(helper_after.encode()),
            whole_helper_inverse_restores_parent=True, zero_private_destination_payload_words=40,
            zero_private_destination_tag_words=8, full8_A_B_words_and_tags_unchanged=True,
            original_unmasked_scope_and_masked_copyback_byte_identical=True,
            original_destination_parameter_retained_unused=True),
        compiler_private_storage=basis['compiler_private_storage'],
        source_proof=[
            'Only each stage sole EvaluateOneStatic declaration gains inline, plus the helper include path. ISPC inline forces expansion according to the exact pinned local guide. All399 layers, generic evaluator, layer/header bounds and9380 ordered tuples retain exact bytes.',
            'Within the existing unmasked snapshot block only private destination slot0 changes to integer zero: all40 payload words and8 tag words. SourceA/B payload and tags retain all physical8 lanes and exact words, including inactive/unknown metadata.',
            'Root and arithmetic_review source plus existing emitted-write review established that all5 unchanged tagged C++ functions compute active destination payload/tags solely from A/B. OriginalD appears only in integer inactive-lane merges; by-value Triple computation does not feed that merge back into write_tag.',
            'Canonical helper arguments remain d0/a1/b2. Private inactiveD is discarded by unchanged original-mask copyback; active outputs/tags are fully written. Source alias snapshots remain before any original-bank write. Original Add/Subtract/Multiply identity shortcuts still precede packing.',
            'Captured mask, existing unmasked scope, C++call placement/gating, full A/B snapshots and masked copyback are unchanged. No mask0 skip, unmasked C++call, FP operation/control write, arithmetic/radius/formula change or original destination publication is added.',
            'Generic fallback, original scalar/v3 paths, exact eight-lane ABI, original AQ symbol mappings and strict target flags remain. Same explicit576B private arrays; no KernelContext/host/shader buffer/count/cache/governor change.'],
        required_postcompile_gate=dict(before_any_link_or_runtime=True,
            expected_removed_EvaluateOneStatic_functions=16,
            reject_remaining_EvaluateOneStatic_symbols_or_calls=True,
            inspect_all_executable_sections=True, reject_unknown_or_non_v4_ISA_FMA_x87_FP_control_writes=True,
            record_actual_code_and_stack_changes=True, require_original_finite_and_adoption_gates=True,
            no_SSA_stack_reuse_or_performance_inference=True),
        scope='Source-only next structural candidate after preserved AR frame failure. Three forced opcode inline qualifiers and integer zero-initialization of private D only. All original semantics/counts/limits and remaining scope barriers stay. Prospective3raw compiles only; no compile/link/provider/workload/timing or qualification transfer. Parent static16scope elimination is the future emitted-code decision; full SSA, stack reuse, spills and speed are unproven.')
    with (W / 'derivation.json').open('x') as stream:
        stream.write(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(source_only=True, raw_sources=3, opcode_inline_qualifiers=3,
        private_D_zero_words=48, literal_instruction_calls=9380, compiled=False, executed=False)))


if __name__ == '__main__':
    assert sys.argv[1:] == ['derive']
    derive()
