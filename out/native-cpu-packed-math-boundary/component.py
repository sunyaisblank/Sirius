"""Derive a private packed boundary; never invoke a compiler or workload."""
from pathlib import Path
import hashlib
import json
import re
import sys
import zipfile

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-packed-math-boundary'
STAGES = ('transport', 'endpoint', 'dopri_phase')
OPERATIONS = (('Add', 'add', 2), ('Multiply', 'multiply', 4),
              ('Subtract', 'subtract', 3), ('Divide', 'divide', None), ('Sqrt', 'sqrt', None))
PINS = {
    'out/native-cpu-v4-raw-packets/derivation.json': 'f50b0e9c73272be4a7945fb88a2c8437cba1c5df8c62e6450354987cd1f2aa8a',
    'out/native-cpu-v4-raw-packets/component.py': 'a9a2b06be656a5f15bab970d3906a5099ce295a9cf0ba2111e85fe47f9ee80b5',
    'out/native-cpu-local-identities/identity_store.isph': '0d464b6ce820abc38301d8f305ba6802ae04e5c13eff80c953bad118b786be2c',
    'out/native-cpu-value-fact-stores/value_facts.isph': '7a245fb084f620a394f7920723b55d66b6c5986a21af78c78ea92cfb149c9f23',
    'out/native-cpu-value-fact-stores/retained8.cpp': '60455d362b1c98b30e9975c35011c6a193ffcf98a5872ddb45036c0e5364961d',
    'out/native-cpu-v4-arithmetic/optional.cpp': 'bcd7804fc7935de78846bf6317d9545b7df78afa934957d74691e74927b028d5',
    'out/native-cpu-v4-raw-provider/provider/source-contract.json': 'edd1786861ae446611b6b90fd1b48d1c7b67a9bfa4318d921d0e85efc664bf3d',
    'out/native-cpu-v4-raw-provider/provider/linked-native-proof.json': '73c9fb917bbfdad6c79b75722448141ab61f5f686db61a5ba2b0953cc1ce2ea4',
    'attestations/qualification/984ba457/native-cpu-v4-raw-packets-finite-compatibility/MANIFEST.json': 'b41bc17d44c4757c51415a1282da22692c85cbd589e92c11724e47da01cdc19b',
    'bin/linux-gcc/native-cpu-simd-rows/compiler/ispc.html': '11ca3d2a0584083270fb77cb2dd95ca6d048ec0f5b3ac25a9638935a2c9652fc',
}
PACKING = '''// Compiler-private marshalling only; the target remains exactly eight lanes.
// Uniform indices retain the original static guards. All source words, including
// inactive physical lanes, were initialized by the unchanged full-context memset.
// The temporary d/a/b slots are 0/1/2. No original bank pointer reaches C++.
static inline void SiriusPackedMathSnapshot_0(
    uniform uint32_t * uniform bank, uniform uint32_t * uniform tags,
    uniform uint32_t destination, uniform uint32_t ai, uniform uint32_t bi,
    uniform uint32_t * uniform packed, uniform uint32_t * uniform packedTags) {
    unmasked {
        for (uniform uint32_t component = 0U; component < 5U; ++component) {
            packed[component*programCount+programIndex] =
                bank[(destination*5U+component)*programCount+programIndex];
            packed[(5U+component)*programCount+programIndex] =
                bank[(ai*5U+component)*programCount+programIndex];
            packed[(10U+component)*programCount+programIndex] =
                bank[(bi*5U+component)*programCount+programIndex];
        }
        packedTags[programIndex] = tags[destination*programCount+programIndex];
        packedTags[programCount+programIndex] = tags[ai*programCount+programIndex];
        packedTags[2U*programCount+programIndex] = tags[bi*programCount+programIndex];
    }
    // ISPC restores the incoming mask here. There is no external call above.
}

static inline void SiriusPackedMathCopyBack_0(
    uniform uint32_t * uniform bank, uniform uint32_t * uniform tags,
    uniform uint32_t destination,
    uniform uint32_t * uniform packed, uniform uint32_t * uniform packedTags) {
    // The incoming execution mask is unchanged; no unmasked stores to the bank.
    for (uniform uint32_t component = 0U; component < 5U; ++component)
        bank[(destination*5U+component)*programCount+programIndex] =
            packed[component*programCount+programIndex];
    tags[destination*programCount+programIndex] = packedTags[programIndex];
}
'''


def doc(path):
    return json.loads(path.read_text())


def data_pin(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def key(path):
    return str(path.relative_to(R)) if path.is_relative_to(R) else str(path)


def function(text, name):
    match = re.search(r'^static noinline void ' + name + r'\([^{}]*\)\n\{', text, re.M)
    assert match, name
    end, depth = match.end(), 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}')
        end += 1
    return text[match.start():end]


def derive():
    assert W == R / 'out/native-cpu-packed-math-boundary'
    assert not (W / 'derivation.json').exists() and not (W / 'packinghelper.isph').exists()
    assert not B.exists()
    for name, sha in PINS.items():
        assert pin(R / name)['sha256'] == sha, name
    basis_path = R / 'out/native-cpu-v4-raw-packets/derivation.json'
    basis = doc(basis_path)
    provider = doc(R / 'out/native-cpu-v4-raw-provider/provider/source-contract.json')
    linked = doc(R / 'out/native-cpu-v4-raw-provider/provider/linked-native-proof.json')
    assert linked['passed'] and len(linked['raw_v4_objects']) == 3
    assert provider['raw_component_manifest'] in PINS
    manifest = doc(R / provider['raw_component_manifest'])
    archive = R / manifest['archive']['path']
    assert pin(archive) == {k: manifest['archive'][k] for k in ('bytes', 'sha256')}
    with zipfile.ZipFile(archive) as zipped:
        archived = zipped.read('out/native-cpu-v4-raw-packets/derivation.json')
        assert data_pin(archived) == pin(basis_path)
        assert data_pin(archived) == manifest['members']['out/native-cpu-v4-raw-packets/derivation.json']
    compiler = R / basis['compiler']['path']
    assert pin(compiler) == {k: basis['compiler'][k] for k in ('bytes', 'sha256')}
    guide = R / 'bin/linux-gcc/native-cpu-simd-rows/compiler/ispc.html'
    guide_text = guide.read_text()
    assert '// and again only the a &lt; b instances' in guide_text
    assert 'not\nmade if none of the program instances need to make the call' in guide_text
    inputs = {name: pin(R / name) for name in PINS}
    inputs[key(archive)] = pin(archive)
    inputs[key(compiler)] = pin(compiler)
    inputs[key(Path(__file__).resolve())] = pin(Path(__file__).resolve())
    helper = W / 'packinghelper.isph'
    include = '#include "' + str(helper) + '"\n'
    sources, stages, commands = {helper: PACKING.encode()}, {}, []
    total = 0
    for slug in STAGES:
        old_path = R / basis['stages'][slug]['derived_path']
        assert pin(old_path) == basis['derived'][key(old_path)]
        source = old_path.read_text()
        with zipfile.ZipFile(archive) as zipped:
            assert zipped.read(key(old_path)) == source.encode()
        assert source.startswith(basis['stages'][slug]['prefix'])
        assert source.count('memset((void * uniform)&context, 0, sizeof(context));') == 1
        assert re.search(r'uint32_t retainedScratch_0\[\d+\];', source)
        assert re.search(r'uint32_t SiriusNormalized_0\[\d+\];', source)
        candidate = source
        changes = []
        for label, operation, identity in OPERATIONS:
            name = 'RT' + label + 'Store_0'
            old = function(candidate, name)
            signature, body = old.split('{', 1)
            target = 'sirius_native_tagged_' + operation + '8'
            old_call = ('    ' + target + '((uniform uint32_t * uniform)&SiriusContext->retainedScratch_0[0],\n'
                '        (uniform uint32_t * uniform)&SiriusContext->SiriusNormalized_0[0],\n'
                '        SiriusDestination, SiriusA, SiriusB, (uniform uint32_t)lanemask());\n')
            assert body.endswith(old_call + '}'), (slug, name)
            first = body[:-len(old_call + '}')]
            if identity is None:
                assert first == '\n'
            else:
                assert first == ('\n    if(SiriusTryIdentity_0((uniform uint32_t * uniform)&SiriusContext->retainedScratch_0[0],\n'
                    '        (uniform uint32_t * uniform)&SiriusContext->SiriusNormalized_0[0],\n'
                    '        SiriusDestination,SiriusA,SiriusB,' + str(identity) + 'U)) return;\n')
            # Only the bounded integer snapshot is unmasked, in the included helper.
            # The original external call remains in the incoming execution mask.
            addition = '''    uniform uint32_t SiriusPackedMask = (uniform uint32_t)lanemask();
    uniform uint32_t SiriusPackedBank[120];
    uniform uint32_t SiriusPackedTags[24];
    SiriusPackedMathSnapshot_0(
        (uniform uint32_t * uniform)&SiriusContext->retainedScratch_0[0],
        (uniform uint32_t * uniform)&SiriusContext->SiriusNormalized_0[0],
        SiriusDestination, SiriusA, SiriusB, SiriusPackedBank, SiriusPackedTags);
'''
            addition += '    ' + target + '(SiriusPackedBank, SiriusPackedTags,\n        0U, 1U, 2U, SiriusPackedMask);\n'
            addition += '''    SiriusPackedMathCopyBack_0(
        (uniform uint32_t * uniform)&SiriusContext->retainedScratch_0[0],
        (uniform uint32_t * uniform)&SiriusContext->SiriusNormalized_0[0],
        SiriusDestination, SiriusPackedBank, SiriusPackedTags);
}'''
            new = signature.replace('static noinline void ', 'static inline void ', 1) + '{' + first + addition
            assert candidate.count(old) == 1
            candidate = candidate.replace(old, new, 1)
            changes.append((old, new))
        anchor = '// Direct private-scratch result ABI; original arithmetic helpers retained.\n'
        assert candidate.count(anchor) == 1
        candidate = candidate.replace(anchor, include + anchor, 1)
        inverse = candidate.replace(include, '', 1)
        for old, new in reversed(changes):
            assert inverse.count(new) == 1
            inverse = inverse.replace(new, old, 1)
        assert inverse == source
        literals = lambda text: [line for line in text.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
        assert literals(candidate) == literals(source)
        assert len(literals(source)) == basis['stages'][slug]['literal_instruction_calls']
        total += len(literals(source))
        destination = B / 'ispc-source' / (slug + '.ispc')
        sources[destination] = candidate.encode()
        command, = [row for row in basis['raw_compilations'] if row['stage'] == slug]
        argv = list(command['argv'])
        assert argv[1] == str(old_path)
        assert argv[2:-3] == ['-O2', '--target=avx512skx-x8', '--cpu=skx', '--arch=x86-64',
                              '--addressing=32', '--PIC', '--opt=disable-fma', '--math-lib=system']
        argv[1], argv[-1] = str(destination), str(B / 'objects' / (slug + '.o'))
        assert argv[-3:-1] == ['--emit-obj', '-o']
        commands.append(dict(stage=slug, argv=argv, output=argv[-1]))
        stages[slug] = dict(basis_path=key(old_path), basis=pin(old_path), derived_path=key(destination),
            literal_instruction_calls=len(literals(source)), whole_source_inverse_restores_AQ=True,
            changed_functions=[dict(name='RT' + label + 'Store_0', original=data_pin(old.encode()),
                candidate=data_pin(new.encode()), identity_shortcut_bytes_unchanged=identity is not None)
                for (label, _op, identity), (old, new) in zip(OPERATIONS, changes)],
            preprocessor_symbol_map=basis['stages'][slug]['preprocessor_symbol_map'],
            EvaluateStaticLayer_and_EvaluateOneStatic_scopes_unchanged=True,
            original_export_and_all_external_declarations_unchanged=True,
            original_compile_argv=command['argv'])
        inputs[key(old_path)] = pin(old_path)
    assert total == 9380
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    derived = {}
    for path, data in sources.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            stream.write(data)
        derived[key(path)] = pin(path)
    reference_raw = [row['object'] for row in linked['raw_v4_objects']]
    reference_paths = {row['path'] for row in reference_raw}
    reused = provider['native_objects'][:7] + [row for row in provider['adapter_input_objects']
                                               if row['path'] not in reference_paths]
    receipt = dict(prepared_only=True, compiled=False, executed=False, source_inputs=inputs,
        derived=derived, stages=stages, raw_compilations=commands,
        literal_instruction_calls_preserved=total, compiler=basis['compiler'],
        raw_component_manifest=provider['raw_component_manifest'],
        reused_objects=reused, reference_optional_raw_objects=reference_raw,
        compiler_private_storage=dict(payload_words=120, tag_words=24, bytes=576,
            slots=dict(destination=0, a=1, b=2), physical_lanes=8,
            stack_reuse_established=False, host_allocations_changed=False, KernelContext_layout_changed=False),
        mask_and_FP_proof=[
            'Original Add/Multiply/Subtract integer identity shortcut executes first with byte-identical source; its return still bypasses packing and C++.',
            'Capture lanemask before the unmasked snapshot. Only uniform-index integer loads/stores occur there; all original physical scratch/tag lanes are initialized by the unchanged full-context memset.',
            'The pinned ISPC guide restores the incoming mask after unmasked. The same extern C call occurs outside unmasked, retains original any-active call gating and receives the captured original mask. No new mask0 skip or forced C call exists.',
            'Destination, both operands and all three tags are snapshotted for all8 lanes before C++ or any original-bank store. This preserves inactive operand words used by full-vector C++ arithmetic and its FP flags.',
            'The unchanged C++ helper reads complete source operands before destination writes; fixed distinct temporary slots preserve original d/a/b aliases and ai==bi values. Its inactive destination merges see exact original destination words/tags.',
            'Only destination payload and destination tag are copied back, under the restored original execution mask. Unknown/legacy metadata, signed zeros, invalid stored results and unchanged source/input words retain their original meanings.',
            'Sqrt still fills temporary B/tagB from original SiriusB; the original canonical call passes B=A. Unchanged C++ tagged_sqrt ignores b_index, so fixed b=2 remains a filled, unused slot.',
            'Marshalling performs no FP operations or FP-control writes; the direct v4 symbol mappings and original strict target flags are unchanged. Exact generated-code/mask/flag equivalence remains a future qualification gate.'],
        local_semantics=dict(guide=dict(path=key(guide), **pin(guide)),
            sections=['re-establishing-the-execution-mask', 'interoperability-with-the-application'],
            source_review_GO=True, compiler_or_runtime_validation=False),
        required_postcompile_gate=dict(basis['required_postcompile_gate'],
            inspect_actual_wrapper_inlining=True, record_stack_and_code_size=True,
            prove_original_banks_do_not_escape_to_external_math=True,
            compare_exact_whole_bank_alias_masks_metadata_and_FP_flags=True,
            original_provider_cohort_benefit_controls_frame_gates_still_required=True),
        scope='Source-only local escape-narrowing hypothesis. Five shared RT*Store bodies become inline; one private integer-packing include is added. All9380 ordered tuples, layer/evaluator scopes, exact8lane public ABI, AQ export/helper symbol mappings, arithmetic source, FP/ISA flags, private KernelContext layout, host/shader buffers, counts, policy, scheduling, cache and governors are unchanged. Each general helper call introduces explicit compiler-private576B temporary storage; stack reuse, SROA/SSA, spills, speed and frame completion are unproven. No compiler, link, provider, workload, timing, compatibility transfer or scientific/release acceptance. Coupling to another separately qualified arithmetic component is a later owner decision.')
    with (W / 'derivation.json').open('x') as stream:
        stream.write(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(source_only=True, raw_sources=3, changed_store_functions=15,
        literal_instruction_calls=total, private_temporary_bytes=576, compiled=False, executed=False)))


if __name__ == '__main__':
    assert sys.argv[1:] == ['derive']
    derive()
