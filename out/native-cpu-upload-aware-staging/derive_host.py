"""Derive upload-aware staging with exact inverse to the matched7f2 source."""
from pathlib import Path
import hashlib
import json
import re
import sys

root = Path.cwd()
work = Path(__file__).resolve().parent
source_path = root / 'out/native-cpu-row-ownership/retained_compute.cpp'
destination = root / 'bin/linux-gcc/native-cpu-upload-aware-staging/source/retained_compute.cpp'
assert sys.argv[1:] in ([], ['--write'])


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        assert sys.argv[1:] == ['--write'], path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)


copy_proof_path = root / 'out/native-cpu-host-copy-prepare/proof.json'
copy_proof = json.loads(copy_proof_path.read_text())
source = source_path.read_text()
assert pin(source.encode()) == copy_proof['sources'][str(source_path.relative_to(root))]['basis']
header_path = root / 'src/sirius/backend/retained_compute.h'
header = header_path.read_text()
assert header.count('bool program_uploaded = false;') == 1
assert header.index('  private:') < header.index('    struct Stage {')
# Both upload sites publish the flag only after the full WriteBuffer succeeds.
assert re.findall(r'program_uploaded\s*=\s*(\w+)', source) == ['true', 'true']
assert source.count('program_uploaded') == 4
dispatch_write = '''    if (stage.program_uploaded) input = input.first(1 + active_rows * stage.input_words_per_row);
    const auto write_started = Clock::now();
    auto status = device_.WriteBuffer(stage.buffers[0], std::as_bytes(input));
    const auto write_finished = Clock::now();
    stage.stats.write_buffer_ms += milliseconds(write_started, write_finished);
    if (!status) return status;
    stage.program_uploaded = true;'''
pair_write = '''        if (stage.program_uploaded) input = input.first(1 + active * stage.input_words_per_row);
        const auto started = Clock::now();
        auto status = device_.WriteBuffer(stage.buffers[0], std::as_bytes(input));
        stage.stats.write_buffer_ms += milliseconds(started, Clock::now());
        if (status) {
            stage.program_uploaded = true;'''
assert source.count(dispatch_write) == source.count(pair_write) == 1
prepare = source[source.index('base::Expected<void> RetainedCompute::PrepareSoftwareRendererStages('):
                 source.index('std::uint64_t RetainedCompute::RequiredBufferBytes(')]
assert 'program_uploaded' not in prepare and '.input' not in prepare and '->input' not in prepare
assert 'auto status = write_header(0);' in prepare
assert 'auto restored = write_header(static_cast<std::uint32_t>(capacity_));' in prepare
reset = source[source.index('void RetainedCompute::ResetStatistics()'):
               source.index('base::Expected<void> RetainedCompute::Dispatch(')]
assert 'program_uploaded' not in reset and '.input' not in reset
# The only other host reader is guarded by the last successful Step's issuance.
owner_guard = '''            if (token.issuer_.lock() != owner.identity || token.generation_ != owner.generation ||
                token.row_ != row || !owner.issued[row])
                continue;
            const auto* original = transport_.input.data() + 1 + row * 230;'''
assert source.count(owner_guard) == 1
assert source.count('std::fill(owner.issued.begin(), owner.issued.end(), false);') == 1
step = source[source.index('base::Expected<std::vector<RetainedStepOutput>> RetainedCompute::Step('):
              source.index('base::Expected<std::vector<RetainedEndpointOutput>> RetainedCompute::Endpoint(')]
assert step.index('InvalidateStageSevenOwner();') < step.index('std::memcpy(')
assert step.count('stage_seven_->issued[row] = result[row].valid && result[row].rhs_valid;') == 1
assert step.index('for (std::size_t row = 0; row < result.size(); ++row)') < step.index('stage_seven_->issued[row]')
# Freeze the complete spelling-level inventory, including constructor and writes.
input_uses = [line.strip() for line in source.splitlines() if re.search(r'(?:\.|->)input\b', line)]
assert len(input_uses) == 24
assert sum('std::fill_n(' in line for line in input_uses) == 9
assert sum('std::memcpy(' in line for line in input_uses) == 9
assert [line for line in input_uses if not any(word in line for word in
    ('std::fill_n(', 'std::memcpy('))] == [
    'stage.input.resize(bytes[0] / 4);',
    'stage.input[0] = static_cast<std::uint32_t>(capacity);',
    'stage.input.begin() + 1 + layout.input_words * capacity);',
    'const auto* original = transport_.input.data() + 1 + row * 230;',
    'auto input = std::span(stage.input);',
    'auto input = std::span(stage.input);']
spacing_patch = root / 'out/native-cpu-host-representation/host.patch'
lines = spacing_patch.read_text().splitlines(True)
assert sum(line.startswith('@@') for line in lines) == 1
hunk = lines[next(i for i, line in enumerate(lines) if line.startswith('@@')) + 1:]
assert all(line[:1] in (' ', '+', '-') for line in hunk)
old_spacing = ''.join(line[1:] for line in hunk if line[:1] in (' ', '-'))
new_spacing = ''.join(line[1:] for line in hunk if line[:1] in (' ', '+'))
assert source.count(old_spacing) == 1
text = source.replace(old_spacing, new_spacing, 1)
fallback = ('return double(std::nextafter(magnitude, std::numeric_limits<float>::infinity())) -\n'
            '               double(magnitude);')
assert old_spacing.count(fallback) == new_spacing.count(fallback) == 1
edits, sites = [], []
pattern = re.compile(r'    std::fill_n\((\w+_)\.input.begin\(\) \+ 1, capacity_ \* (\d+), 0U\);\n'
                     r'    std::memcpy\(\1\.input.data\(\) \+ 1, (\w+)\.data\(\), \3\.size_bytes\(\)\);')
for match in pattern.finditer(text):
    stage, stride, span = match.groups()
    replacement = (f'    if (!{stage}.program_uploaded)\n'
                   f'        std::fill_n({stage}.input.begin() + 1 + {span}.size() * {stride},\n'
                   f'                    (capacity_ - {span}.size()) * {stride}, 0U);\n'
                   f'    std::memcpy({stage}.input.data() + 1, {span}.data(), {span}.size_bytes());')
    edits.append((match[0], replacement))
    sites.append({'stage': stage, 'words_per_row': int(stride), 'input': span})
assert len(edits) == 9
assert sites == [{key: value for key, value in item.items() if key != 'line'}
                 for item in copy_proof['sources'][str(source_path.relative_to(root))]['sites']]
for old, new in edits:
    assert text.count(old) == 1
    text = text.replace(old, new, 1)
restored = text
for old, new in reversed(edits):
    assert restored.count(new) == 1
    restored = restored.replace(new, old, 1)
assert restored.count(new_spacing) == 1
restored = restored.replace(new_spacing, old_spacing, 1)
assert restored == source
inputs = {str(path.relative_to(root)): pin(path.read_bytes())
          for path in (source_path, header_path, spacing_patch, copy_proof_path,
                       root / 'out/native-cpu-host-copy-prepare/matched-host.patch')}
libm = Path('/usr/lib/x86_64-linux-gnu/libm.so.6')
inputs[str(libm)] = pin(libm.read_bytes())
proof = {'basis_path': str(source_path.relative_to(root)), 'basis': pin(source.encode()),
    'derived_path': str(destination.relative_to(root)), 'derived': pin(text.encode()),
    'source_inputs': inputs, 'full_inverse_restores_matched_source': True,
    'normal_spacing_patch_applied_once_in_private_copy': True,
    'original_spacing_fallback_verbatim': True, 'inactive_only_input_zero_ranges': sites,
    'inactive_clear_guards': 9, 'inactive_clear_only_before_first_successful_upload': True,
    'all_other_source_bytes_unchanged': True, 'source_selectors_and_all_shader_arrays_unchanged': True,
    'normal_spacing_scope': 'Original abs remains first. Only finite positive normal magnitudes below maxfinite use the exact exponent-derived normal binary64 ULP. All zero/subnormal/maxfinite/nonfinite cases retain the original nextafter expression and its exception/errno behavior. Root independently reviews/tests linked-libm FEnv/errno equivalence before runtime admission.',
    'copy_scope': 'Complete active memcpy retained at all nine sites. Inactive host rows are cleared only while that Stage.program_uploaded is false. Header/program host suffix, every consumed upload byte, device state, authorized host reads and refusal/rollback behavior are unchanged. Inactive unconsumed host bytes may differ after a successful full upload.',
    'upload_wire_invariant': {
        'flag_default_false': True, 'flag_assignment_count': 2,
        'assignments_only_true_after_successful_write': True,
        'stage_input_reference_lines': input_uses,
        'false': 'Full serialization follows unchanged inactive clear plus full active memcpy. Failed writes leave the flag false, so a smaller retry clears newly inactive rows before its next full write; any partial-write behavior receives the same wire bytes and call sequence.',
        'true': 'Successful full upload already installed the immutable program tail. Every later serialization is exactly header plus active rows; each active row is fully copied before upload, including newly active rows after growth. Skipped inactive host bytes never enter those uploads. A later dispatch/read/decode/timing failure does not invalidate the completed upload.',
        'paired_endpoint_dense': 'The two Stage flags remain independent. Endpoint write success followed by Dense write failure leaves Endpoint true and Dense false, exactly as before; write order, early returns and each subsequent serialized range are unchanged.',
        'preparation_and_reset': 'Preparation writes only device header zero/restored capacity and never changes the host input or uploaded-program flag. ResetStatistics changes neither. Cancellation and errors retain original ordering and flag lifetime.',
        'stage_seven_host_reader': 'Step invalidates every issued bit before validation/staging. Only successfully decoded active result rows may be issued; authorization checks matching issuer/generation/row and issued before reading that transport input row. No inactive stale host tail is consumed.',
        'all_nine_active_copies_and_all_other_source_bytes_retained': True,
        'no_new_fields_allocations_upload_ranges_or_dispatches': True},
    'limits': 'Derived source only; no compiler/build/runtime, bottleneck, speedup, scientific or release claim.'}
publish(destination, text.encode())
publish(work / 'host-source-proof.json', (json.dumps(proof, indent=2) + '\n').encode())
assert source_path.read_text() == source
print(json.dumps({'private_source_derived': True, 'zero_range_changes': 9,
                  'spacing_changes': 1, 'full_inverse_exact': True, 'compiled_or_executed': False}))
