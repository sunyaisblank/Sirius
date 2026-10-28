"""Derive source only: a dispatch capability skips repeated program comparisons."""
from pathlib import Path
import difflib
import hashlib
import json

R = Path(__file__).resolve().parents[2]
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-dispatch-admission-token'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
PINS = {
    "attestations/qualification/e17097b9/native-cpu-inline-private-bank-finite-compatibility/MANIFEST.json": {
        "bytes": 47027,
        "sha256": "757366c154706de8085e8c5c00be7d0fc99dfca8faefd8905905c5d2b7ec74a8"
    },
    "attestations/qualification/e17097b9/native-cpu-inline-private-bank-finite-compatibility/component.zip": {
        "bytes": 10462956,
        "sha256": "bbe723c238042fe256722d82a1032d838da7eb9a96cb3afbcc1e530e72eadb64"
    },
    "bin/linux-gcc/native-cpu-endpoint-header-domain/admission.o": {
        "bytes": 194776,
        "sha256": "0a9a768f2605e70c7be24ed3c8e3177b2e7301dfced4b7eb002d6a38fc8d3d20"
    },
    "bin/linux-gcc/native-cpu-endpoint-header-domain/endpoint.o": {
        "bytes": 263152,
        "sha256": "7730ce37e8eb50ec8e83abd9430401e6faf99e2ad4bc8bf1f4bda4c36b486847"
    },
    "bin/linux-gcc/native-cpu-exact-product-zero-residual/optional.o": {
        "bytes": 256280,
        "sha256": "31594a31204d85e8fec6ef1463807c5c1117eddc6fe4e8ad297f4f02362ed59f"
    },
    "bin/linux-gcc/native-cpu-inline-private-bank/ispc-source/dopri_phase.ispc": {
        "bytes": 403852,
        "sha256": "f9c1197a5c06814705d6bf7f46cf0e8ce217d7f42534627fadbfa208f51dac8b"
    },
    "bin/linux-gcc/native-cpu-inline-private-bank/ispc-source/endpoint.ispc": {
        "bytes": 795302,
        "sha256": "d8ae7cdde50211537678ec197f18c374edd2cec653974d9d0b7c33c50db8c1e1"
    },
    "bin/linux-gcc/native-cpu-inline-private-bank/ispc-source/transport.ispc": {
        "bytes": 626976,
        "sha256": "0fc91ff0b756fddd055d5c0322f5eece0557a64799ef31893933be33cdba10f6"
    },
    "bin/linux-gcc/native-cpu-inline-private-bank/objects/dopri_phase.o": {
        "bytes": 2438712,
        "sha256": "5f6a79dd19817031ac6d180c46761f317f4bc0f53b908865fd528a8b768c5ab1"
    },
    "bin/linux-gcc/native-cpu-inline-private-bank/objects/endpoint.o": {
        "bytes": 5018192,
        "sha256": "39a33f821780caf92fa435da33e23cea659fa59f62dcdd3e0087821461e88a90"
    },
    "bin/linux-gcc/native-cpu-inline-private-bank/objects/transport.o": {
        "bytes": 3732232,
        "sha256": "3cbd66c1e54076f7c66552f5b42f504229a2a6a39bf9969973fba5d3e89c2d06"
    },
    "bin/linux-gcc/native-cpu-inline-private-provider/rows_adapter.o": {
        "bytes": 11657104,
        "sha256": "d505a833db55e2a3435a9ce5e38598f947e731b9d640f48771048739eccd5be6"
    },
    "bin/linux-gcc/native-cpu-local-identities/objects/dopri_phase.o": {
        "bytes": 1193360,
        "sha256": "0dc138c31434c91f1af671f81513a9a5ce6e35a20eab05846a7c12e5b19e766d"
    },
    "bin/linux-gcc/native-cpu-local-identities/objects/endpoint.o": {
        "bytes": 934344,
        "sha256": "1330cf6875d0db24aaad70f407da8607285d651cbc38e6d471c8247e6e398a54"
    },
    "bin/linux-gcc/native-cpu-local-identities/objects/transport.o": {
        "bytes": 1042176,
        "sha256": "04bc4c45bf21b7df05ac4fb95283aba9e7b208cd73d336b83093bb26c33b8728"
    },
    "bin/linux-gcc/native-cpu-scalar-zero/objects/dopri_phase.o": {
        "bytes": 138152,
        "sha256": "ffb626bc3898cd6f860f077378f68e4c800f2f772fbbb878fec3eff580dbeb95"
    },
    "bin/linux-gcc/native-cpu-scalar-zero/objects/transport.o": {
        "bytes": 206512,
        "sha256": "06d6769d09ed9bbabc8fce7b97748358f8721bbdcab8fa4ec15dba6718184f81"
    },
    "bin/linux-gcc/native-cpu-scalar-zero/typed-batch/dopri_phase.cpp": {
        "bytes": 2866,
        "sha256": "e481ae98846a7d8d509e2e3e17e8b8d72ea561e30929a165c85e9026c57d93d6"
    },
    "bin/linux-gcc/native-cpu-scalar-zero/typed-batch/dopri_phase.native.hpp": {
        "bytes": 410475,
        "sha256": "0808c37704c7f98423bb1b5b0f57c74feef5d27f015b74caf4dc62c6255ca626"
    },
    "bin/linux-gcc/native-cpu-scalar-zero/typed-batch/endpoint.native.hpp": {
        "bytes": 813886,
        "sha256": "b95b630d7fdb7a50db8f3551d3e48a3c932b9d5944dabc4b45afe301cc45c698"
    },
    "bin/linux-gcc/native-cpu-scalar-zero/typed-batch/transport.cpp": {
        "bytes": 2925,
        "sha256": "fe8daf2bfc2d4b0e87370b6b559c70f3687166667609c8fadcc54cbe69452083"
    },
    "bin/linux-gcc/native-cpu-scalar-zero/typed-batch/transport.native.hpp": {
        "bytes": 635340,
        "sha256": "e9d605e77065dd22c7018526c4d8d7bec27960d5c8d2be27306d22f0545d4ea4"
    },
    "bin/linux-gcc/native-cpu-slang-bridge/typed-batch/dopri_phase.program.hpp": {
        "bytes": 44967,
        "sha256": "3df9a067a4e65fd58959fd84076bb6c743e33fb15f26aa898d699029c6cd64fe"
    },
    "bin/linux-gcc/native-cpu-slang-bridge/typed-batch/endpoint.program.hpp": {
        "bytes": 85168,
        "sha256": "b6f8ba278b7e9df65860143d9711972603b4f8a1850b2ade138bf17f625150d5"
    },
    "bin/linux-gcc/native-cpu-slang-bridge/typed-batch/transport.program.hpp": {
        "bytes": 64160,
        "sha256": "9c581d3a98843f9cd831a7814e8b7ad77a0f40ab2c8cb29301df67cc6bad1ffc"
    },
    "bin/linux-gcc/native-cpu-v4-arithmetic/dispatch.o": {
        "bytes": 4704,
        "sha256": "d7b9aa5ae74d4b791b9e6df59d0e2c10d3ae5fd079aeacf15b80a899e40b6411"
    },
    "bin/linux-gcc/native-cpu-v4-arithmetic/v3.o": {
        "bytes": 238552,
        "sha256": "0af01ed41b0189bc047b10b7d853a0e4ebdbb5c6eb02f729870401c2eb24f287"
    },
    "bin/linux-gcc/native-cpu-v4-raw-packets/objects/rows_adapter.o": {
        "bytes": 4760,
        "sha256": "c9ee8c3c860673d0db9014b8d08b5cf95a38d4b999a6e74ce9da409f3c832d62"
    },
    "out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8/src/gallium/drivers/llvmpipe/lp_native_cpu_bridge.h": {
        "bytes": 6421,
        "sha256": "ffa014c8d06e66c33935477fa11a13b775b9cc479860a3be5d6943adc7c3568f"
    },
    "out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8/src/gallium/drivers/llvmpipe/lp_state_cs.c": {
        "bytes": 109640,
        "sha256": "34b1f093c09543c7f1c48873081d014849088b3839316dc6b3671d2b6add9cd1"
    },
    "out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8/src/gallium/drivers/llvmpipe/sirius_native_bridge.h": {
        "bytes": 1939,
        "sha256": "a6f3da671bd9507f9f97300266ac0feb0d34a179a659143151d4e7e33d60f971"
    },
    "out/native-cpu-active-prefix/native_packet.hpp": {
        "bytes": 2445,
        "sha256": "9df6d65d1c89016694a7f242a8fdd61979db6826c4f0d2619b2b8954a9ca6eee"
    },
    "out/native-cpu-endpoint-header-domain/admission.cpp": {
        "bytes": 2165,
        "sha256": "5dd4535e222d8e82546627ea0174e3f31b3deb45648790da92ade40f5499585d"
    },
    "out/native-cpu-endpoint-header-domain/build-plan.json": {
        "bytes": 2146,
        "sha256": "c2b824c6046d3ede2eb9d0984a2f9c7978c0fa6052184915ff24ba1793aac373"
    },
    "out/native-cpu-endpoint-header-domain/endpoint.cpp": {
        "bytes": 2915,
        "sha256": "5d54fbc79103a44bf017eb30e8a2be9e80e1ff93557b89c6622b83ade10bf025"
    },
    "out/native-cpu-endpoint-header-domain/native_packet.hpp": {
        "bytes": 3767,
        "sha256": "9da2464ac8812752e7091ccad7e53ff7c276e2666d6a9750f18aa4c66bab91b7"
    },
    "out/native-cpu-inline-private-provider/provider/source-contract.json": {
        "bytes": 127856,
        "sha256": "96d8ef11e0ae93f341cf18498e93947dd3b99a3878852581f432aac9b66a3b15"
    },
    "out/native-cpu-scalar-zero/cpp-compile-result.json": {
        "bytes": 14560,
        "sha256": "0fd97b1d7d8680b3e63358f827e52f28528f88b65f11a3d265adc905ef01f1db"
    },
    "out/native-cpu-simd-packet-policy/rows_bridge.h": {
        "bytes": 1538,
        "sha256": "d52b7cb2080ee35f0f06ea0df9e5e71f8694a3cc45024b486b5cb7a50b9dc282"
    },
    "out/native-cpu-slang-bridge/native_bridge.h": {
        "bytes": 1055,
        "sha256": "57d78d5d8a61cdd04f4d9e3301d600c96fc716c3581c7d2cdfb045956798b243"
    },
    "out/native-cpu-v4-raw-packets/derivation.json": {
        "bytes": 129718,
        "sha256": "f50b0e9c73272be4a7945fb88a2c8437cba1c5df8c62e6450354987cd1f2aa8a"
    },
    "out/native-cpu-v4-raw-packets/rows_adapter.cpp": {
        "bytes": 4831,
        "sha256": "8452185cf367983695af42404dbc0415fc5f8b497fe5a56da317aa6a8e7bd99f"
    }
}
STAGES = ('transport', 'endpoint', 'dopri_phase')
LAYOUTS = {'transport': (230, 4699, 15967), 'endpoint': (225, 3829, 21143),
           'dopri_phase': (1810, 4149, 10689)}
SCALARS = {
    'transport': 'bin/linux-gcc/native-cpu-scalar-zero/typed-batch/transport.cpp',
    'endpoint': 'out/native-cpu-endpoint-header-domain/endpoint.cpp',
    'dopri_phase': 'bin/linux-gcc/native-cpu-scalar-zero/typed-batch/dopri_phase.cpp',
}


def pin(path):
    data = path.read_bytes()
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def doc(path):
    return json.loads(path.read_text())


def replace_one(text, before, after):
    assert text.count(before) == 1, before
    return text.replace(before, after, 1)


def extract_function(text, name):
    start = text.index('extern "C" bool ' + name + '(')
    opened = text.index('{', start)
    depth = 1
    end = opened + 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}')
        end += 1
    return text[start:end]


def derive():
    assert PINS
    for name, expected in PINS.items():
        assert pin(R / name) == expected, name
    sources = {}
    header_basis = (R / 'out/native-cpu-endpoint-header-domain/native_packet.hpp').read_text()
    old_template = 'template<bool endpoint_header_domain = false>'
    # Endpoint's default-false specialization must retain the original strict
    # Transport/DP guard too, despite sharing the newer Endpoint header file.
    strict_compare = '    return std::memcmp(suffix, program, program_words * sizeof(std::uint32_t)) == 0;'
    begin = header_basis.index('    if constexpr (endpoint_header_domain) {')
    end = header_basis.index(strict_compare, begin)
    strict_projection = header_basis[:begin] + header_basis[end:]
    strict_projection = replace_one(strict_projection, '#include <endpoint.program.hpp>\n', '')
    strict_projection = replace_one(strict_projection, old_template + '\n', '')
    assert strict_projection == (R / 'out/native-cpu-active-prefix/native_packet.hpp').read_text()
    new_template = ('// dispatch_program_admitted is private: the identical live views and immutable\n'
                    '// program have passed complete admission in this synchronous Mesa job.\n'
                    '// It skips only program-byte comparisons; every other guard below remains.\n'
                    'template<bool endpoint_header_domain = false, bool dispatch_program_admitted = false>')
    header = replace_one(header_basis, old_template, new_template)
    endpoint_compare = '        return std::memcmp(suffix, program, register_word * sizeof(std::uint32_t)) == 0 &&'
    skip_endpoint = '        if constexpr (dispatch_program_admitted) return true;\n'
    skip_strict = '    if constexpr (dispatch_program_admitted) return true;\n'
    header = replace_one(header, endpoint_compare, skip_endpoint + endpoint_compare)
    header = replace_one(header, strict_compare, skip_strict + strict_compare)
    inverse = replace_one(header, skip_endpoint + endpoint_compare, endpoint_compare)
    inverse = replace_one(inverse, skip_strict + strict_compare, strict_compare)
    inverse = replace_one(inverse, new_template, old_template)
    assert inverse == header_basis
    sources['native_packet.hpp'] = header

    bridge = '''#pragma once

// Private entry contract: only lp_native_cpu_packet_mode's positive job-local
// native_preflight_passed branch may call rows_admitted_packet. That preflight
// validated the same stage, live views/extents, immutable capacity/grid and all
// required canonical program bytes before any output write. The input and
// descriptor table remain immutable through the existing synchronous pool wait.
// Output is disjoint from the whole input and descriptor table. Each worker
// still checks all original pointer/extent/header/shape/alias/packet bounds.
// Arithmetic refusals remain handled execution, never late fallback.
// Public try_packet/admit_packet entries remain fully guarded; no token survives
// this job and no pointer/global/TLS cache or caller-supplied bool is introduced.
#include "../native-cpu-simd-packet-policy/rows_bridge.h"

#ifdef __cplusplus
extern "C" {
#define SIRIUS_ADMITTED_NOEXCEPT noexcept
#else
#define SIRIUS_ADMITTED_NOEXCEPT
#endif
#define SIRIUS_ADMITTED_ARGUMENTS \\
    const void *input, uint32_t input_bytes, void *output, uint32_t output_bytes, \\
    uint32_t first_x, uint32_t y, uint32_t z, uint32_t grid_x, uint32_t grid_y, \\
    uint32_t grid_z, uint32_t groups
'''
    for stage in STAGES:
        for suffix in ('admitted_packet', 'rows_admitted_packet'):
            bridge += ('__attribute__((visibility("hidden")))\n'
                       f'bool sirius_native_{stage}_{suffix}(SIRIUS_ADMITTED_ARGUMENTS) SIRIUS_ADMITTED_NOEXCEPT;\n')
    bridge += '''#undef SIRIUS_ADMITTED_ARGUMENTS
#undef SIRIUS_ADMITTED_NOEXCEPT
#ifdef __cplusplus
}
#endif
'''
    sources['admitted_bridge.h'] = bridge
    scalar_proofs = {}
    for stage in STAGES:
        original = (R / SCALARS[stage]).read_text()
        public = 'sirius_native_' + stage + '_try_packet'
        private = 'sirius_native_' + stage + '_admitted_packet'
        old = extract_function(original, public)
        new = replace_one(old, 'extern "C" bool ' + public + '(',
                          'extern "C" __attribute__((visibility("hidden"))) bool ' + private + '(')
        before = 'sirius_cpu_bridge::admit_packet' + ('<true>' if stage == 'endpoint' else '') + '('
        after = 'sirius_cpu_bridge::admit_packet<' + ('true' if stage == 'endpoint' else 'false') + ', true>('
        new = replace_one(new, before, after)
        restored = replace_one(new, after, before)
        restored = replace_one(restored,
            'extern "C" __attribute__((visibility("hidden"))) bool ' + private + '(',
            'extern "C" bool ' + public + '(')
        assert restored == old
        extra = ('\n#include "admitted_bridge.h"\n\n' + new + '\n\n'
                 f'static_assert(std::is_same_v<decltype(&{private}), decltype(&{public})>);\n')
        sources[stage + '.cpp'] = original + extra
        assert sources[stage + '.cpp'][:-len(extra)] == original
        scalar_proofs[stage] = dict(public_wrapper_exact_source_prefix=True,
            private_body_inverse_restores_public_function=True,
            shared_original_generated_body_and_types=True, guard=after[:-1])

    adapter_basis = (R / 'out/native-cpu-v4-raw-packets/rows_adapter.cpp').read_text()
    macro_start = adapter_basis.index('#define ROWS_ARGUMENTS')
    macro_end = adapter_basis.index('// ISPC\'s input pointer type', macro_start)
    macros = adapter_basis[macro_start:macro_end]
    extra = ('\n#include "admitted_bridge.h"\n#include "native_packet.hpp"\n'
             '#include <transport.program.hpp>\n#include <dopri_phase.program.hpp>\n\n' + macros)
    adapter_proofs = {}
    for stage in STAGES:
        public = 'sirius_native_' + stage + '_rows_try_packet'
        private = 'sirius_native_' + stage + '_rows_admitted_packet'
        old = extract_function(adapter_basis, public)
        new = replace_one(old, 'extern "C" bool ' + public + '(',
                          'extern "C" __attribute__((visibility("hidden"))) bool ' + private + '(')
        new = replace_one(new, 'sirius_native_' + stage + '_try_packet(ROWS_VALUES)',
                          'sirius_native_' + stage + '_admitted_packet(ROWS_VALUES)')
        before = 'sirius_native_' + stage + '_admit_packet(ROWS_VALUES)'
        iw, ow, words = LAYOUTS[stage]
        after = (f'sirius_cpu_bridge::admit_packet<{str(stage == "endpoint").lower()}, true>(ROWS_VALUES, '
                 f'{iw}U, {ow}U,\n        sirius_cpu_bridge::{stage}_program, {words}U)')
        new = replace_one(new, before, after)
        restored = replace_one(new, after, before)
        restored = replace_one(restored, 'sirius_native_' + stage + '_admitted_packet(ROWS_VALUES)',
                               'sirius_native_' + stage + '_try_packet(ROWS_VALUES)')
        restored = replace_one(restored,
            'extern "C" __attribute__((visibility("hidden"))) bool ' + private + '(',
            'extern "C" bool ' + public + '(')
        assert restored == old
        extra += '\n' + new + '\n'
        extra += f'static_assert(std::is_same_v<decltype(&{private}), decltype(&{public})>);\n'
        adapter_proofs[stage] = dict(private_inverse_restores_public_function=True,
            full8_v3_v4_ISA_and_scalar_tail_selection_unchanged=True)
    extra += '\n#undef ROWS_VALUES\n#undef ROWS_ARGUMENTS\n'
    sources['rows_adapter.cpp'] = adapter_basis + extra
    assert sources['rows_adapter.cpp'][:-len(extra)] == adapter_basis

    name = 'src/gallium/drivers/llvmpipe/lp_native_cpu_bridge.h'
    original = (S / name).read_text()
    anchor = '   switch (job->native_cpu_stage) {\n   case 1:\n      return sirius_native_transport_rows_try_packet'
    branch = ('   /* Complete preflight belongs to this exact live synchronous job.\n'
              '    * Only the repeated immutable-program comparison is omitted below;\n'
              '    * private entries retain every original view/header/packet guard.\n'
              '    * Existing ordinary/refusal/JIT paths remain fully guarded.\n'
              '    */\n'
              '   if (!preflight && job->native_preflight_passed) {\n'
              '      switch (job->native_cpu_stage) {\n')
    for index, stage in enumerate(STAGES, 1):
        branch += (f'      case {index}:\n'
                   f'         return sirius_native_{stage}_rows_admitted_packet(input->u, input->num_elements,\n'
                   '            (void *)output->u, output->num_elements, x, y, z,\n'
                   '            job->grid_size[0], job->grid_size[1], job->grid_size[2], groups);\n')
    branch += '      }\n   }\n'
    patched = replace_one(original, anchor, branch + anchor)
    assert replace_one(patched, branch + anchor, anchor) == original
    changes = {name: patched}
    name = 'src/gallium/drivers/llvmpipe/sirius_native_bridge.h'
    original = (S / name).read_text()
    anchor = '// Pure runtime ISA/OS-state query; scalar baseline remains supported.\n'
    declarations = ('// Private dispatch-admitted workers: only the existing job-local positive\n'
                    '// preflight branch may call these; all public entries stay fully guarded.\n')
    for stage in STAGES:
        declarations += (f'__attribute__((visibility("hidden"))) bool sirius_native_{stage}_rows_admitted_packet'
                         '(SIRIUS_NATIVE_PACKET_ARGUMENTS) SIRIUS_NATIVE_NOEXCEPT;\n')
    declarations += '\n'
    patched = replace_one(original, anchor, declarations + anchor)
    assert replace_one(patched, declarations + anchor, anchor) == original
    changes[name] = patched
    sources['mesa.patch'] = ''.join(''.join(difflib.unified_diff(
        (S / name).read_text().splitlines(True), candidate.splitlines(True),
        fromfile='a/' + name, tofile='b/' + name)) for name, candidate in changes.items())

    old_scalar_commands = doc(R / 'out/native-cpu-scalar-zero/cpp-compile-result.json')['objects']
    endpoint_command = doc(R / 'out/native-cpu-endpoint-header-domain/build-plan.json')['commands'][1]
    commands = []
    for stage in STAGES:
        args = list(endpoint_command if stage == 'endpoint' else
                    next(row['argv'] for row in old_scalar_commands if row['stage'] == stage))
        args[args.index('-c') + 1] = str(W / (stage + '.cpp'))
        args[args.index('-o') + 1] = str(B / (stage + '.o'))
        commands.append(dict(role=stage, argv=args))
    args = list(doc(R / 'out/native-cpu-v4-raw-packets/derivation.json')['adapter_compile_argv'])
    args[args.index('-c') + 1] = str(W / 'rows_adapter.cpp')
    args[args.index('-o') + 1] = str(B / 'rows_adapter.o')
    commands.append(dict(role='rows_adapter', argv=args))
    contract = doc(R / 'out/native-cpu-inline-private-provider/provider/source-contract.json')
    reused = contract['native_objects'][3:7] + contract['adapter_input_objects'][1:]
    for row in reused:
        assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}, row['path']
    for name, expected in PINS.items():
        assert pin(R / name) == expected, name
    for name, text in sources.items():
        with (W / name).open('x') as stream:
            stream.write(text)
    proof = dict(prepared_only=True, compiled=False, executed=False,
        source_inputs=dict(PINS, **{str(Path(__file__).relative_to(R)): pin(Path(__file__))}),
        derived={name: pin(W / name) for name in sources}, prospective_compile_commands=commands,
        new_compilations_if_selected=4, reused_objects=reused,
        scalar_source_inverse=scalar_proofs, adapter_source_inverse=adapter_proofs,
        guard_inverse_restores_AP_header=True, default_false_projection_restores_original_strict_guard=True,
        Mesa_patch_inverse_restores_exact_two_headers=True,
        existing_public_wrapper_sources_preserved=True, original_pure_admission_object_reused=True,
        no_job_layout_change=True, no_new_state_buffers_allocations_or_caches=True,
        raw3_generated_sources_objects_and_all_ISPC_static_generic_guards_unchanged=True,
        capability=dict(authority='lp_cs_job_info.native_preflight_passed',
            establishment='memset-false; positive pure preflight at x=0/groups=1 before scheduling or output writes',
            binding='same job stage/current descriptor views, full extents, capacity/header, grid and canonical suffix; no caller-supplied bool',
            lifetime='existing direct task or lp_cs_tpool_wait_for_task before job stack expires',
            immutable='whole input incl capacity/program and descriptor table remain unchanged during the synchronous dispatch',
            nonalias='original guards require output disjoint from whole input and descriptor table; active row outputs occupy disjoint original strides',
            repeated_guards='all original pointer/alignment/overflow/extent/capacity/grid/y/z/groups/first_x checks retained in private entries',
            FP_ISA='original cs_exec_fn save/IEEE32/restore, v3+OS query, v4+OS selector, scalar/tail/generated row code unchanged',
            failure='failed preflight never sets capability and builds original JIT before writes; false-worker abort invariant unchanged',
            zero_header='unchanged separate zero-active no-op schedules no workers and grants no capability'),
        endpoint_program=dict(words=21143, compared_words=21142, exempt_relative_word=11241,
            allowed_value='any uint32 unchanged; selected specialized row refusals remain generated-shader decisions'),
        guard_calls=dict(current_SIMD='1+floor(G/8)+G%8 program comparisons per successful native dispatch',
            current_scalar='1+G', candidate='one pure preflight comparison; every worker retains non-program guards'),
        source_findings='Only repeated whole-program comparisons are elided on the private positively admitted path. No cross-dispatch pointer cache. Input/descriptors are read-only during the current synchronous job; output cannot overwrite them. All preflight word refusals happen before output publication; scalar and SIMD workers retain all other original guards. Exact shader-internal program selection/static/generic guards remain untouched. Source inverse/contract proof only; no emitted code, link, execution, timing, provider, cohort, scientific or release acceptance transferred.')
    with (W / 'derivation.json').open('x') as stream:
        json.dump(proof, stream, indent=2)
        stream.write('\n')
    print(json.dumps(dict(source_only=True, derived_files=len(sources), prospective_compilations=4,
                         Mesa_patch_files=len(changes), source_inverses_passed=True)))


if __name__ == '__main__':
    derive()
