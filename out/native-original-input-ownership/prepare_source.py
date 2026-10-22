"""Prepare one private complete-input ownership candidate; no compiler or dispatch."""
from pathlib import Path
import ast, difflib, hashlib, json, shutil
R = Path.cwd(); OLD = R / 'out/native-endpoint-phase-sharing'
W = Path(__file__).resolve().parent; B = R / 'bin/linux-gcc/native-original-input-ownership'
def replace(text, old, new, count=1):
    assert text.count(old) == count, (old[:100], text.count(old), count)
    return text.replace(old, new)
def pin(p): return {'bytes': p.stat().st_size, 'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
def dump(p, d): p.write_text(json.dumps(d, indent=2) + '\n')
assert not (W / 'compiler-result.json').exists() and not (W / 'actions.json').exists()
for name in ['shaders', 'ordinary', 'readers', 'runtime']: (W / name).mkdir(exist_ok=True)
for name in ['modules', 'matched']: (B / name).mkdir(parents=True, exist_ok=True)
for name in ['owned_guard.py', 'run_owned.py', 'build_matched.py', 'make_build_actions.py']:
    text = (OLD / name).read_text().replace('native-endpoint-phase-sharing', W.name)
    if name == 'build_matched.py':
        text = text.replace('Two existing native32 selectors only', 'Three existing native32 selectors only')
        # Freeze this generated selector before compilation, as in the prior candidate.
    if name == 'make_build_actions.py':
        text = text.replace('exactly one changed native32 Endpoint module', 'exactly two changed native32 Transport/Endpoint modules')
        text = text.replace('private two-module selector', 'private three-module selector')
    (W / name).write_text(text)
for name in ['providers.json', 'provider-reuse.json']:
    shutil.copyfile(OLD / name, W / name)
for name in ['provider-history', 'reused']:
    shutil.copytree(OLD / name, W / name, dirs_exist_ok=True)
observer = OLD.parent.parent / 'bin/linux-gcc/native-endpoint-phase-sharing/matched/audited-observer.o'
shutil.copyfile(observer, B / 'matched/audited-observer.o')
assert pin(observer) == pin(B / 'matched/audited-observer.o')

# Only native32 receives these private sources. Keep every other array untouched.
transport = (R / 'src/sirius/kernels/retained_transport.slang').read_text()
transport = replace(transport, 'RetainedTriple ReadOriginal(uint index, uint row) {', '''#if defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_FP64) || defined(SIRIUS_RETAINED_FMA32) || defined(SIRIUS_RETAINED_STAGE7_OWNER)
#error Complete input ownership belongs only to the existing native binary32 route.
#endif
groupshared uint originalTerms[46 * 5];
RetainedTriple ReadOriginalInput(uint index, uint row) {''')
mark = 'RetainedTriple ReadInput(uint index, uint row, uint rowBase) {'
transport = replace(transport, mark, '''void CacheOriginal(uint index, uint row) {
    RetainedTriple r = ReadOriginalInput(index, row);
    uint base = index * 5;
    originalTerms[base] = RPBits(r.hi);
    originalTerms[base + 1] = RPBits(r.lo);
    originalTerms[base + 2] = RPBits(r.tail);
    originalTerms[base + 3] = RPBits(r.error);
    originalTerms[base + 4] = r.valid;
    if (r.valid == 0) RejectRetainedProgram();
}
RetainedTriple ReadOriginal(uint index, uint row) {
    uint base = index * 5;
    RetainedTriple r;
    r.hi = RPFromBits(originalTerms[base]);
    r.lo = RPFromBits(originalTerms[base + 1]);
    r.tail = RPFromBits(originalTerms[base + 2]);
    r.error = RPFromBits(originalTerms[base + 3]);
    r.valid = originalTerms[base + 4];
    return r;
}
''' + mark)
transport = replace(transport, '    if (rowWords < 2404) return;', '    if (rowWords < 2404) return;\n    if (lane == 0) retainedStatus = 0;')
transport = replace(transport, '''    for (uint i = 0; i < 46; ++i)
        if (ReadOriginal(i, row).valid == 0) return;''', '''    // Shape gates and the initial output/status publication precede packing.
    // Each immutable field has one owner; publish all five words before use.
    for (uint i = lane; i < 46; i += SIRIUS_RETAINED_EXECUTION_LANES)
        CacheOriginal(i, row);
    GroupMemoryBarrierWithGroupSync();
    if (retainedStatus != 0) return;''')
(W / 'shaders/retained_transport.slang').write_text(transport)
endpoint = (OLD / 'shaders/retained_endpoint.slang').read_text()
endpoint = replace(endpoint, '#if defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_FP64)', '#if defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_FP64) || defined(SIRIUS_RETAINED_FMA32) || defined(SIRIUS_RETAINED_STAGE7_OWNER)')
endpoint = replace(endpoint, 'groupshared uint nativeProjection[7];', '''groupshared uint nativeProjection[7];
groupshared uint originalTerms[45 * 5];
groupshared uint firstInvalidOriginal;''')
entry_prefix, entry_body = endpoint.rsplit('#else\n[shader("compute")]', 1)
entry_body = replace(entry_body, '    if (rowWords < 769) return;', '    if (rowWords < 769) return;\n    if (lane == 0) firstInvalidOriginal = 45;')
endpoint = entry_prefix + '#else\n[shader("compute")]' + entry_body
old_loop = '''    for (uint i = 0; i < 45; ++i) {
        uint address = 1 + 225 * row + 5 * i;
        RetainedTriple value =
            RTPack(RPFromBits(inputs[address]), RPFromBits(inputs[address + 1]),
                   RPFromBits(inputs[address + 2]), RPFromBits(inputs[address + 3]));
        if (inputs[address + 4] != 1 || value.valid == 0) return;
        if (lane == 0) WriteStored(rowBase + 504 + 5 * i, value);
    }
#ifndef SIRIUS_RETAINED_PORTABLE
    AllMemoryBarrierWithGroupSync();
#endif'''
endpoint = replace(endpoint, old_loop, '''    // The initial output join publishes sentinel45. Normalize each immutable
    // original once, but commit only the original ordered prefix on refusal.
    for (uint i = lane; i < 45; i += SIRIUS_RETAINED_EXECUTION_LANES) {
        uint address = 1 + 225 * row + 5 * i, base = 5 * i;
        RetainedTriple value =
            RTPack(RPFromBits(inputs[address]), RPFromBits(inputs[address + 1]),
                   RPFromBits(inputs[address + 2]), RPFromBits(inputs[address + 3]));
        if (inputs[address + 4] != 1) value.valid = 0;
        originalTerms[base] = RPBits(value.hi);
        originalTerms[base + 1] = RPBits(value.lo);
        originalTerms[base + 2] = RPBits(value.tail);
        originalTerms[base + 3] = RPBits(value.error);
        originalTerms[base + 4] = value.valid;
        if (value.valid == 0) {
            uint previous;
            InterlockedMin(firstInvalidOriginal, i, previous);
        }
    }
    GroupMemoryBarrierWithGroupSync();
    for (uint i = lane; i < firstInvalidOriginal; i += SIRIUS_RETAINED_EXECUTION_LANES) {
        uint base = 5 * i;
        RetainedTriple value;
        value.hi = RPFromBits(originalTerms[base]);
        value.lo = RPFromBits(originalTerms[base + 1]);
        value.tail = RPFromBits(originalTerms[base + 2]);
        value.error = RPFromBits(originalTerms[base + 3]);
        value.valid = originalTerms[base + 4];
        WriteStored(rowBase + 504 + base, value);
    }
    AllMemoryBarrierWithGroupSync();
    if (firstInvalidOriginal != 45) return;''')
(W / 'shaders/retained_endpoint.slang').write_text(endpoint)

# Private builder: include complete caches, dedicated status, and the already
# accepted seven-word projection in the conservative16KiB source bound.
builder = (R / 'scripts/build-retained-kernels.py').read_text()
builder = replace(builder, 'stage_seven_owner=False):', 'stage_seven_owner=False, complete_inputs=False):')
builder = replace(builder, '    original_inputs = ({', '''    complete_input_words = 0
    if complete_inputs:
        if fp64 or portable or fma or stage_seven_owner or source.stem not in ("retained_transport", "retained_endpoint"):
            raise ValueError("complete input ownership is native32 Transport/Endpoint only")
        complete_input_words = 230 if source.stem == "retained_transport" else 225 + 1
        if source.stem == "retained_endpoint":
            projection_words = 7
    original_inputs = ({''')
builder = replace(builder, 'coefficients * 5 + projection_words) * 4 + 8', 'coefficients * 5 + projection_words + complete_input_words) * 4 + 8')
(W / 'builder.py').write_text(builder)

modules = (OLD / 'build_modules.py').read_text().replace('native-endpoint-phase-sharing', W.name)
modules = replace(modules, '"""Build only the two ownership candidates with the original shader authority."""', '"""Build two complete-input candidates; reuse the admitted DopriPhase module."""')
modules = replace(modules, "R / 'scripts/build-retained-kernels.py'", "W / 'builder.py'")
start = modules.index("        if Path(argv[1]).stem == 'retained_endpoint':")
end = modules.index('    commands.append(argv)', start)
modules = modules[:start] + '''        assert Path(argv[1]).stem in ['retained_endpoint', 'retained_transport']
        argv[1] = str(W / 'shaders' / Path(argv[1]).name)
        argv[2:2] = ['-I', str(R / 'src/sirius/kernels')]
''' + modules[end:]
modules = replace(modules, "    ('Endpoint', 'endpoint', programs.build_endpoint_program),", "    ('Transport', 'transport', programs.build_transport_program),\n    ('Endpoint', 'endpoint', programs.build_endpoint_program),")
modules = replace(modules, '    specialized = None', '    specialized = None; transport_specialized = None')
modules = replace(modules, "    text = (C / 'retained_kernels.h').read_text()", '''    if kind == 'Transport':
        sp = programs.build_schwarzschild_transport_program(parallel=True)
        specialized_words = [sp['instructions'], sp['registers'], len(sp['outputs'])]
        specialized_words += sp['outputs'] + sp['operations'] + sp['layer_offsets']
        transport_specialized = (len(encoded), len(specialized_words), len(sp['layer_offsets']) - 1)
        encoded += specialized_words
    text = (C / 'retained_kernels.h').read_text()''')
modules = replace(modules, '        endpoint_specialized=specialized)', '        endpoint_specialized=specialized, transport_specialized=transport_specialized, complete_inputs=True)')
modules = replace(modules, '    (B / (\'owned-\' + stem + \'.spvasm\')).write_text(assembly)', '''    # Count the actually emitted workgroup variables, including every scalar.
    constants = {k: int(v) for k, v in re.findall(r'(%\\S+) = OpConstant %uint (\\d+)', assembly)}
    types = {k: 4 for k in re.findall(r'(%\\S+) = OpTypeInt 32 [01]', assembly)}
    for key, element, length in re.findall(r'(%\\S+) = OpTypeArray (%\\S+) (%\\S+)', assembly):
        if element in types and length in constants: types[key] = types[element] * constants[length]
    pointers = dict(re.findall(r'(%\\S+) = OpTypePointer Workgroup (%\\S+)', assembly))
    shared = {key: types[pointers[pointer]] for key, pointer in re.findall(r'(%\\S+) = OpVariable (%\\S+) Workgroup', assembly)}
    assert sum(shared.values()) <= 16384
    expected = (plan['registers'] * 5 + (125 + 230 + 1 if kind == 'Transport' else 225 + 1 + 7 + 1)) * 4
    assert sum(shared.values()) == expected, (shared, expected)
    (B / ('owned-' + stem + '.spvasm')).write_text(assembly)''')
modules = replace(modules, "        'scalar_algorithm_unchanged': True}", "        'scalar_algorithm_unchanged': True, 'actual_shared_variables_bytes': shared, 'actual_shared_bytes': sum(shared.values())}")
modules = replace(modules, "'compiled_modules': 1", "'compiled_modules': 2")
(W / 'build_modules.py').write_text(modules)

for name in ['capture.py', 'capture_features.py', 'capture_jit.py', 'feature-layout.json', 'float-layout.json',
             'owned_guard.py', 'properties_capture.py', 'read_features.py', 'run_owned.py', 'make_actions.py']:
    text = (OLD / 'ordinary' / name).read_text().replace('native-endpoint-phase-sharing', W.name)
    if name == 'make_actions.py':
        text = text.replace("[('endpoint','Endpoint')]", "[('transport','Transport'),('endpoint','Endpoint')]")
        text = text.replace("inputs.add('bin/linux-gcc/" + W.name + "/modules/owned-endpoint.spv')", "inputs |= {'bin/linux-gcc/" + W.name + "/modules/owned-endpoint.spv', 'bin/linux-gcc/" + W.name + "/modules/owned-transport.spv'}")
    (W / 'ordinary' / name).write_text(text)
for name in ['native_evidence.py', 'numeric_gates.py', 'decimal_check.py', 'expected-readback-layout.json', 'run_abba.py', 'read_results.py', 'read_controls.py']:
    text = (OLD / 'readers' / name).read_text().replace('native-endpoint-phase-sharing', W.name)
    if name == 'read_results.py':
        text = text.replace("native['endpoint']='kEndpointNativeOwnedShader';", "native['transport']='kTransportNativeOwnedShader';native['endpoint']='kEndpointNativeOwnedShader';")
        text = text.replace("{'endpoint':'Endpoint','dopri_phase':'DopriPhase'}", "{'transport':'Transport','endpoint':'Endpoint','dopri_phase':'DopriPhase'}")
        text = text.replace("stem in ['endpoint','dopri_phase']", "stem in ['transport','endpoint','dopri_phase']")
        text = text.replace('under two ownership modules', 'under three ownership modules')
    if name == 'read_controls.py':
        text = text.replace("native.update(endpoint=", "native.update(transport='kTransportNativeOwnedShader', endpoint=")
        text = text.replace("for kind in ['Endpoint', 'DopriPhase']:", "for kind in ['Transport', 'Endpoint', 'DopriPhase']:")
        start = text.index("        if key == '/usr/lib/x86_64-linux-gnu/liblzma.so.5.4.5':")
        end = text.index('        assert pin(p)', start)
        text = text[:start] + text[end:]
        text = text.replace("       'postexecution_library_identity_restoration': doc(W.parent / 'library-restoration/identity-restoration.json'),\n", '')
    (W / 'readers' / name).write_text(text)
for name in ['prepare.py', 'run_case.py', 'run_owned.py', 'owned_guard.py']:
    text = (OLD / 'runtime' / name).read_text().replace('native-endpoint-phase-sharing', W.name)
    if name == 'prepare.py':
        text = text.replace("'bin/linux-gcc/" + W.name + "/modules/owned-endpoint.spv',", "'bin/linux-gcc/" + W.name + "/modules/owned-endpoint.spv','bin/linux-gcc/" + W.name + "/modules/owned-transport.spv',")
        text = text.replace('same two-module ownership candidate', 'same three-module ownership candidate')
        text = text.replace('native32 Endpoint owns one complete curved projection and DopriPhase distributes all362 original checks', 'native32 Transport owns all46 normalized inputs, Endpoint preserves the exact ordered valid prefix through its first-invalid join and its single complete projection, and DopriPhase distributes all362 original checks')
    (W / 'runtime' / name).write_text(text)
source = (R / 'src/sirius/backend/retained_compute.cpp').read_text()
for kind in ['Transport', 'Endpoint', 'DopriPhase']:
    source = replace(source, 'shader(k' + kind + 'Shader, k' + kind + 'Fp64Shader',
                     'shader(k' + kind + 'NativeOwnedShader, k' + kind + 'Fp64Shader')
(W / 'retained_compute.cpp').write_text(source)
for p in W.rglob('*.py'): ast.parse(p.read_text(), filename=str(p))
patches = []
for old, new in [(R / 'src/sirius/kernels/retained_transport.slang', W / 'shaders/retained_transport.slang'),
                 (OLD / 'shaders/retained_endpoint.slang', W / 'shaders/retained_endpoint.slang'),
                 (R / 'scripts/build-retained-kernels.py', W / 'builder.py')]:
    patches.extend(difflib.unified_diff(old.read_text().splitlines(True), new.read_text().splitlines(True), str(old.relative_to(R)), str(new.relative_to(R))))
(W / 'prototype.patch').write_text(''.join(patches))
dump(W / 'source-scope.json', {'basis_source': 'd9dbd65796aab5dd36e17481198ae1f6b6201a2c',
    'preservation': '1bdfb8a8e8bb5a0865e1913c7febf483db48fe35',
    'change': 'Private native32 complete immutable input ownership in Transport/Endpoint; keep accepted single projection and Dopri module.',
    'Transport_cache_previously_rejected': 'attestations/native-vulkan/3045067/transport-immutable-input-cache-rejected/prototype.patch',
    'material_new_context': 'Actual genuine IEEE32 Mesa/LLVM21 CPU path, current phase-sharing source has all8/all22 PASS, failed cold frame captures Transport/Endpoint dominant submit/wait. Endpoint ordered-prefix ownership is a new mechanism; no speedup is presumed.',
    'Transport_current_phase_fields4_through43_unchanged': True,
    'Endpoint_partial_prefix_preserved_before_first_invalid': True,
    'all_shapes_tables_arithmetic_modes_buffers_workgroups_governors_assertions_preserved': True,
    'observer_reuse': {'source': str(observer.relative_to(R)), 'destination': str((B/'matched/audited-observer.o').relative_to(R)), **pin(observer)},
    'screen_order': 'Two module compiles; actual ordinary emitted/driver code before C8; sole original ABBA22; remaining7 only if22 PASS; original cold frame only ifall8 PASS.',
    'qualification': False})
print(json.dumps({'private_sources_prepared': True, 'new_compilers_or_dispatches': 0, 'public_source_modified': False}))
