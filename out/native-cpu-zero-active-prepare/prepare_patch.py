"""Emit an unapplied two-file zero-header no-op patch and its source proof."""
from pathlib import Path
import difflib
import hashlib
import json
import re

root = Path.cwd()
work = Path(__file__).resolve().parent
mesa = root / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
driver = 'src/gallium/drivers/llvmpipe/'


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        path.write_bytes(data)


def replace_once(text, old, new):
    assert text.count(old) == 1, old[:120]
    return text.replace(old, new, 1)


names = ('lp_native_cpu_bridge.h', 'lp_state_cs.c')
basis_hashes = ('6ff7eb9067b7c418d9111ce3671c6ce0e60239d579c9aa504226d8fd9919b13e',
                '36ac0c4e8b050947656be65cb649390dd23cea2d75c17d1dbe7cb45fd92c8aa7')
basis = {}
for name, expected in zip(names, basis_hashes):
    data = (mesa / driver / name).read_bytes()
    assert pin(data)['sha256'] == expected, name
    basis[name] = data.decode()

bridge = basis['lp_native_cpu_bridge.h']
start = bridge.index('static bool\nlp_native_cpu_packet_mode(')
end = bridge.index('   /* SSBO descriptor extents are bytes.', start)
original_prefix = bridge[start:end]
opening = original_prefix.index('{\n') + 2
signature, extraction = original_prefix[:opening], original_prefix[opening:]
assert extraction.endswith('      return false;\n')
shared_prefix = '''struct lp_native_cpu_views {
   struct lp_jit_buffer input;
   struct lp_jit_buffer output;
};

static bool
lp_native_cpu_extract_views(const struct lp_cs_job_info *job,
                            uint32_t x, uint32_t y, uint32_t z,
                            struct lp_native_cpu_views *views)
{
''' + extraction + '''   views->input = *input;
   views->output = *output;
   return true;
}

''' + signature + '''   struct lp_native_cpu_views views;
   if (!lp_native_cpu_extract_views(job, x, y, z, &views))
      return false;
   const struct lp_jit_buffer *input = &views.input;
   const struct lp_jit_buffer *output = &views.output;
'''
bridge = replace_once(bridge, original_prefix, shared_prefix)
zero_helper = '''
/* Separate from positive native admission: each exact owned SPIR-V returns
 * on header rows==0 before any output/program/FP access. Preserve the same
 * conservative job, descriptor-table, pointer, extent and nonalias domain.
 * No worker, FP scope or positive native packet is needed for this no-op.
 */
static bool
lp_native_cpu_zero_active(const struct lp_cs_job_info *job)
{
   struct lp_native_cpu_views views;
   if (!lp_native_cpu_extract_views(job, 0, 0, 0, &views))
      return false;
   const struct lp_jit_buffer *input = &views.input;
   const struct lp_jit_buffer *output = &views.output;
   const uintptr_t in = (uintptr_t)input->u;
   const uintptr_t out = (uintptr_t)output->u;
   if (!input->u || !output->u ||
       input->num_elements < sizeof(uint32_t) ||
       output->num_elements < sizeof(uint32_t) ||
       ((input->num_elements | output->num_elements) & 3u) ||
       ((in | out) & 3u) ||
       in > UINTPTR_MAX - input->num_elements ||
       out > UINTPTR_MAX - output->num_elements ||
       (in < out + output->num_elements && out < in + input->num_elements))
      return false;
   uint32_t rows;
   memcpy(&rows, input->u, sizeof(rows));
   return rows == 0;
}
'''
bridge += zero_helper
assert bridge.replace(zero_helper, '', 1).replace(shared_prefix, original_prefix, 1) == basis[names[0]]

launch = basis['lp_state_cs.c']
launch_edits = [
    ('   bool native_direct = false;\n',
     '   bool native_direct = false;\n   bool native_zero_active = false;\n'),
    ('''      if (lp_native_cpu_preflight(&job_info)) {
         job_info.native_preflight_passed = true;
         job_info.native_simd_packets = sirius_native_simd_supported();
''', '''      native_zero_active = lp_native_cpu_zero_active(&job_info);
      if (native_zero_active || lp_native_cpu_preflight(&job_info)) {
         /* A zero-header no-op is not positive native packet admission. */
         job_info.native_preflight_passed = !native_zero_active;
         job_info.native_simd_packets = !native_zero_active && sirius_native_simd_supported();
'''),
    ('''   if (job_info.native_simd_packets) {
      scheduled_tasks = job_info.grid_size[0] / 8 + job_info.grid_size[0] % 8;
''', '''   if (native_zero_active) {
      /* num_tasks remains the original logical workgroup count for queries. */
      scheduled_tasks = 0;
   } else if (job_info.native_simd_packets) {
      scheduled_tasks = job_info.grid_size[0] / 8 + job_info.grid_size[0] % 8;
''')]
for old, new in launch_edits:
    launch = replace_once(launch, old, new)
restored = launch
for old, new in reversed(launch_edits):
    restored = replace_once(restored, new, old)
assert restored == basis[names[1]]

# Source proof uses the retained original compiler output; no compiler,
# disassembler, interpreter, shader execution, or provider is launched here.
authority_path = root / 'out/native-cpu-row-ownership/compiler-result.json'
authority = json.loads(authority_path.read_text())
whitelist_path = mesa / 'src/gallium/frontends/lavapipe/sirius_native_modules.h'
whitelist = whitelist_path.read_text()
pipeline_path = mesa / 'src/gallium/frontends/lavapipe/lvp_pipeline.c'
pipeline = pipeline_path.read_text()
assert 'size == sizeof(sirius_native_transport_spirv)' in pipeline
assert 'size == sizeof(sirius_native_endpoint_spirv)' in pipeline
assert 'size == sizeof(sirius_native_dopri_phase_spirv)' in pipeline
source_inputs = {str(p.relative_to(root)): pin(p.read_bytes())
                 for p in (authority_path, whitelist_path, pipeline_path,
                           root / 'out/native-cpu-active-prefix/native_packet.hpp')}
stage_proofs = {}
for stage, slug, stride, minimum, flag, labels in [
    ('Transport', 'transport', 230, 43, '%640',
     ['%712', '%711', '%717', '%716', '%722', '%721', '%730', '%729', '%736']),
    ('Endpoint', 'endpoint', 225, 103, '%flat',
     ['%661', '%660', '%666', '%665', '%671', '%670', '%679', '%678', '%685']),
    ('DopriPhase', 'dopri_phase', 1810, 203, '%445',
     ['%476', '%475', '%481', '%480', '%486', '%485', '%494', '%493', '%500'])]:
    item = authority['modules'][stage]
    module_path = root / item['module']
    module = module_path.read_bytes()
    assert pin(module) == {key: item[key] for key in ('bytes', 'sha256')}
    match, = re.finditer(r'sirius_native_' + slug + r'_spirv\[\] = \{([^}]+)\};', whitelist)
    embedded = bytes(int(value, 16) for value in re.findall(r'0x([0-9a-f]{2})', match[1]))
    assert embedded == module
    shader_path = root / ('out/native-cpu-row-ownership/shaders/retained_' + slug + '.slang')
    shader = shader_path.read_text()
    owner = shader[shader.index('void ComputeOwnedRow(uint3 id, uint lane) {'):]
    marker = ('    if (rows == 0 || rows > 65536 || row >= rows || ni < 1 + ' +
              str(stride) + ' * rows + ' + str(minimum) + ' || no % rows != 0)\n        return;')
    assert marker in owner
    guard_end = owner.index(marker) + len(marker)
    prefix = owner[:guard_end]
    expected_prefix = '''void ComputeOwnedRow(uint3 id, uint lane) {
    uint ni, si, no, so;
    inputs.GetDimensions(ni, si);
    outputs.GetDimensions(no, so);
    if (ni < 1 || id.y != 0 || id.z != 0) return;
    uint rows = inputs[0], row = id.x;
''' + marker
    assert prefix == expected_prefix
    assert 'if (lane == 0) ComputeOwnedRow(id, 0);' in shader
    assembly_path = module_path.with_suffix('.spvasm')
    assembly = assembly_path.read_text()
    owned = assembly[assembly.index('%ComputeOwnedRow = OpFunction'):]
    # The original disassembly makes the rows==0 short circuit explicit.
    # Its true successor propagates a function-private boolean to OpReturn.
    blocks = dict(re.findall(r'(?m)^(%\w+) = OpLabel\n(.*?)(?=^%\w+ = OpLabel\n|\Z)', owned,
                             flags=re.M | re.S))
    # Original labels are padded; normalize only leading text whitespace.
    normalized = '\n'.join(line.strip() for line in owned.splitlines()) + '\n'
    blocks = dict(re.findall(r'(?m)^(%\w+) = OpLabel\n(.*?)(?=^%\w+ = OpLabel\n|\Z)', normalized,
                             flags=re.M | re.S))
    for index, label in enumerate(labels):
        block = blocks[label]
        if index == len(labels) - 1:
            assert block.strip() == 'OpReturn'
        elif index % 2 == 0:
            assert block.strip() == 'OpStore ' + flag + ' %true\nOpBranch ' + labels[index + 1]
        else:
            condition, = re.findall(r'(%\w+) = OpLoad %bool ' + re.escape(flag), block)
            assert re.search(r'OpBranchConditional ' + re.escape(condition) + ' ' +
                             re.escape(labels[index + 1]) + r' %\w+', block)
        assert not any(op in block for op in ('OpFunctionCall', 'OpUMod', 'OpUDiv',
                                              'OpAccessChain', 'OpFAdd', 'OpFMul'))
    first_zero_branch = re.search(r'(%\w+) = OpIEqual %bool (%\w+) %uint_0\n'
                                 r'OpSelectionMerge %\w+ None\nOpBranchConditional \1 ' +
                                 re.escape(labels[0]) + r' %\w+', normalized)
    assert first_zero_branch
    header_pointer, = re.findall(r'(%\w+) = OpAccessChain %_ptr_StorageBuffer_uint %inputs %int_0 %int_0', normalized)
    assert first_zero_branch[2] + ' = OpLoad %uint ' + header_pointer in normalized
    stage_proofs[stage] = {'module': {'path': str(module_path.relative_to(root)), **pin(module)},
        'original_shader': str(shader_path.relative_to(root)), 'zero_header_return_prefix': prefix,
        'original_spirv_assembly': str(assembly_path.relative_to(root)),
        'zero_branch_return_labels': labels,
        'complete_whitelist_module_equal': True,
        'other_63_logical_invocations_return_at_entry': True,
        'zero_path_has_no_output_program_fp_or_division_access': True}
    for p in (module_path, shader_path, assembly_path):
        source_inputs[str(p.relative_to(root))] = pin(p.read_bytes())

new_sources = dict(zip(names, (bridge, launch)))
patch = ''
changed = {}
for name in names:
    relative = driver + name
    old, new = basis[name], new_sources[name]
    patch += ''.join(difflib.unified_diff(old.splitlines(True), new.splitlines(True), n=3,
                                         fromfile='a/' + relative, tofile='b/' + relative))
    changed[relative] = {'basis': pin(old.encode()), 'candidate': pin(new.encode())}
    source_inputs[str((mesa / relative).relative_to(root))] = pin(old.encode())
publish(work / 'provider.patch', patch.encode())
proof = {'patch': pin(patch.encode()), 'changed_files': changed, 'source_inputs': source_inputs,
    'stage_proofs': stage_proofs, 'positive_B_guard_unchanged': True,
    'original_descriptor_extraction_body_shared_byte_exact': True,
    'launch_inverse_restores_basis': True,
    'coarse_native_dispatch_domain_unchanged': True,
    'zero_header_classifier_separate_from_positive_admission': True,
    'input_header_immutable_for_dispatch': True,
    'zero_job_variant': None, 'later_fallback_variant_dirty': 'LP_CSNEW_CS',
    'zero_scheduled_tasks': 0, 'original_logical_query_count': 'uint64_t(num_tasks)*64',
    'worker_fp_scope_or_math_entered_for_zero': False,
    'other_modules_views_and_ineligible_dispatches_keep_original_fallback': True,
    'provider_identity': 'Root must give the integrated provider its next private UUID/cache identity.',
    'scope': 'Unapplied source-only incremental patch. No compiler/build/runtime, numerical/performance/cold-frame/release claim.'}
publish(work / 'source-proof.json', (json.dumps(proof, indent=2) + '\n').encode())
assert all((mesa / driver / name).read_text() == basis[name] for name in names)
print(json.dumps({'source_patch_written': True, 'changed_files': len(changed),
                  'patch_applied': False, 'compiled_or_executed': False}))
