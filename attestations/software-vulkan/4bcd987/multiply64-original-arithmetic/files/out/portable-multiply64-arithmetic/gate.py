"""Original arithmetic gate: root-operated producers and separate terminal reader."""
from pathlib import Path
import collections
import hashlib
import json
import re
import shlex
import struct
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
sys.dont_write_bytecode = True
W = Path(__file__).resolve().parent
R = W.parents[1]
T = R / 'bin/linux-gcc/portable-multiply64-arithmetic'
P = json.loads((W / 'packet.json').read_text())

def record(path):
    path = Path(path)
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': path.stat().st_size, 'sha256': digest,
            'resolved_path': str(path.resolve(strict=True))}

def content(path):
    return {k: v for k, v in record(path).items() if k != 'resolved_path'}

def write(path, value):
    path = Path(path)
    assert not path.exists(), path
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')

def verify_static():
    for path, pin in json.loads((W / 'source-freeze.json').read_text())['pins'].items():
        assert record(path) == pin, path

def corpus_checks():
    fixture = (W / 'portable_binary32_reference.bin').read_bytes()
    assert hashlib.sha256(fixture).hexdigest() == P['original_fixture_sha256']
    assert struct.unpack_from('<III', fixture) == (0x32334250, 1, 174268)
    rows = list(struct.iter_unpack('<6I', fixture[12:]))
    assert len(rows) == 174268
    counts = {str(key): value for key, value in sorted(collections.Counter(row[0] for row in rows).items())}
    assert counts == P['corpus_operation_cardinalities'] and set(counts) == {str(n) for n in range(14)}
    assert (W / 'input.bin').read_bytes() == struct.pack('<I', 174268) + b''.join(
        struct.pack('<4I', op, a, b, 0) for op, a, b, hi, lo, valid in rows)
    assert (W / 'expected.bin').read_bytes() == b''.join(
        struct.pack('<4I', hi, lo, valid, 0) for op, a, b, hi, lo, valid in rows)
    divide = [(a, b) for op, a, b, hi, lo, valid in rows if op == 3]
    finite = [(a, b) for a, b in divide if
              0 < (a & 0x7fffffff) < 0x7f800000 and 0 < (b & 0x7fffffff) < 0x7f800000]
    def significand(word):
        word &= 0x7fffffff
        fraction = word & 0x7fffff
        return (fraction | (1 << 23)) if word >> 23 else fraction << (24 - fraction.bit_length())
    coverage = {'divide_records': len(divide), 'optional_finite_nonzero_branch_records': len(finite),
                'normalization_doubling_records': sum(significand(a) < significand(b) for a, b in finite),
                'normalization_kept_records': sum(significand(a) >= significand(b) for a, b in finite),
                'non_power_of_two_denominator_significand_records': sum(significand(b) != (1 << 23) for a, b in finite)}
    assert coverage == P['division_branch_coverage']
    return {'original_records': 174268, 'operation_cardinalities': counts,
            'division_branch_coverage': coverage,
            'scope': 'Exact original ordered fixture/provider byte joins and source-derived branch partitions; not exhaustive machine proof.'}

def accepted_folder(folder):
    folder = Path(folder)
    owner = json.loads((folder / 'owner.json').read_text())
    assert owner['status'] == 'terminal' and owner['passed'] is True and owner['child_exit'] == 0
    assert owner['stop_reason'] is None and owner['cleanup_errors'] == []
    assert owner['remaining_owned_processes'] == []
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
    return owner

def accepted(label):
    return accepted_folder(W / label)

def relocated(path):
    return Path(P['historical_host_archive']) / 'files' / Path(path).relative_to(R)

def verify_host_reuse():
    history = Path(P['historical_host_root'])
    accepted_folder(history / 'compile-device')
    seal = json.loads((history / 'compile-device/output-seal.json').read_text())
    for key, relative in [('source_freeze', 'source-freeze.json'),
                          ('operator_context', 'operator-context.json'),
                          ('inputs_before', 'compile-device/inputs-before.json')]:
        assert content(history / relative) == {k: seal[key][k] for k in ('bytes', 'sha256')}
    for original, pin in seal['outputs'].items():
        assert content(relocated(original)) == {k: pin[k] for k in ('bytes', 'sha256')}
    assert str(P['device_object']) in [str(relocated(name)) for name in seal['outputs']]
    pair = json.loads(Path(P['historical_pair_review']).read_text())
    assert pair['passed'] and pair['complete_55_enabled_fields_queue_extensions_order_and_feature_chain_equal']
    features = []
    for name in ('stock', 'candidate'):
        folder = history / ('capture-' + name + '-endpoint-normal')
        accepted_folder(folder)
        features.append(json.loads((folder / 'feature-review.json').read_text()))
        assert features[-1]['actual_public_create_arguments_enable_shaderInt64']
        assert features[-1]['successful_logical_device_used_for_actual_complete_shader']
    assert features[0] == features[1] == pair['normalized_feature_report']
    assert len(features[0]['features']) == 55 and features[0]['features']['shaderInt64'] == 1
    assert features[0]['features']['shaderFloat64'] == 1
    return {'device_object': record(P['device_object']), 'original_output_seal': record(history / 'compile-device/output-seal.json'),
            'historical_pair_review': record(P['historical_pair_review']),
            'scope': 'Exact accepted compiled Int64 device object/source and historical public-argument witness reused. No fresh numeric-process argument capture is claimed.'}

def verify_output_seal(label):
    accepted(label)
    seal = json.loads((W / label / 'output-seal.json').read_text())
    assert seal['label'] == label and seal['recipe'] == P['recipes'][label]
    assert seal['source_freeze'] == record(W / 'source-freeze.json')
    assert seal['operator_context'] == record(W / 'operator-context.json')
    assert seal['inputs_before'] == record(W / label / 'inputs-before.json')
    assert seal['outputs'] == {name: record(name) for name in P['outputs'][label]}
    return record(W / label / 'output-seal.json')

def verify_archive(path, expected):
    path = Path(path)
    assert content(path / 'manifest.json') == expected
    manifest = json.loads((path / 'manifest.json').read_text())
    for name, pin in manifest['files'].items():
        item = path / name
        assert not item.is_symlink() and content(item) == pin, item
    return record(path / 'manifest.json')

def verify_reference_acceptances():
    compiler = Path(P['compiler_work']); jit = Path(P['jit_work'])
    compiler_manifest = verify_archive(P['compiler_archive'], P['compiler_archive_manifest'])
    jit_manifest = verify_archive(P['jit_archive'], P['jit_archive_manifest'])
    approval = json.loads((compiler / 'compiler-acceptance.json').read_text())
    assert approval['passed'] is True
    capture = json.loads((jit / 'capture-acceptance.json').read_text())
    terminal = json.loads((jit / 'independent-terminal-review.json').read_text())
    cost = json.loads((jit / 'independent-cost-review.json').read_text())
    assert capture['passed'] and terminal['passed']
    assert capture['complete_55_feature_queue_extension_order_and_chain_equal']
    assert capture['fresh_inferior_and_called_JIT_maps_join_all_three_provider_files']
    assert cost['decision'] == 'COST_PREMISE_SURVIVES_UNPROVEN'
    assert (W / 'candidate/src/sirius/kernels/portable_binary32.h').read_bytes() == (compiler / 'candidate/src/sirius/kernels/portable_binary32.h').read_bytes()
    assert (W / 'stock/src/sirius/kernels/portable_binary32.h').read_bytes() == (compiler / 'stock/src/sirius/kernels/portable_binary32.h').read_bytes()
    return {'compiler_manifest': compiler_manifest, 'compiler_acceptance': record(compiler / 'compiler-acceptance.json'),
            'compiler_result': record(compiler / 'compiler-result.json'), 'jit_manifest': jit_manifest,
            'jit_capture_acceptance': record(jit / 'capture-acceptance.json'), 'jit_terminal_review': record(jit / 'independent-terminal-review.json'),
            'ordinary_cost_review': record(jit / 'independent-cost-review.json'), 'normal_candidate_module': record(P['normal_candidate_module']),
            'scope': 'Immutable compiler and ordinary JIT/capability source receipts only; the short-input canary did not execute multiplication and proves no arithmetic or benefit.'}

def verify_stock_reuse():
    archive = Path(P['stock_archive']); old = Path(P['stock_work']); products = Path(P['stock_build'])
    manifest = verify_archive(archive, P['stock_archive_manifest'])
    frozen = json.loads((old / 'source-freeze.json').read_text())['pins']
    assert len(frozen) == len(P['stock_static_join_rows']) == 737
    for row in P['stock_static_join_rows']:
        assert row['pin'] == frozen[row['original']]
        assert content(row['actual']) == {k: row['pin'][k] for k in ('bytes', 'sha256')}
    result = json.loads((old / 'arithmetic-result.json').read_text())
    independent = json.loads((old / 'independent-terminal-review.json').read_text())
    assert result['passed'] and independent['accepted']
    assert result['original_records'] == 174268 and result['complete_output_words_per_provider'] == 697072
    assert result['host_explicit64_and_both_full_gpu_outputs_equal_original_expected']
    assert len(result['owners']) == 17
    for name, pin in result['owners'].items():
        accepted_folder(old / name)
        assert content(old / name / 'owner.json') == {k: pin[k] for k in ('bytes', 'sha256')}
    owner = accepted_folder(old / 'execute-gpu-stock')
    provider_mappings(owner)
    before = json.loads((old / 'execute-gpu-stock/inputs-before.json').read_text())
    for original, actual in ((R / 'bin/linux-gcc/portable-divide64-arithmetic/gpu', products / 'gpu'),
                             (R / 'bin/linux-gcc/portable-divide64-arithmetic/stock.spv', products / 'stock.spv'),
                             (R / 'out/portable-divide64-arithmetic/input.bin', old / 'input.bin'),
                             (R / 'out/portable-divide64-arithmetic/expected.bin', old / 'expected.bin')):
        key = str(original.relative_to(R))
        assert {k: before[key][k] for k in ('bytes', 'sha256')} == content(actual)
    for name in ('portable_binary32_reference.bin', 'input.bin', 'expected.bin', 'probe.slang', 'gpu.cpp'):
        assert (W / name).read_bytes() == (old / name).read_bytes()
    assert (old / 'execute-gpu-stock/actual.bin').read_bytes() == (W / 'expected.bin').read_bytes()
    rows = [json.loads(line) for line in (old / 'execute-gpu-stock/stdout.log').read_text().splitlines() if line.startswith('{')]
    assert len(rows) == 1 and rows[0]['records'] == 174268 and rows[0]['complete_words'] == 697072 and rows[0]['mismatches'] == 0
    assert rows[0]['driver'] == 'llvmpipe'
    return {'archive_manifest': manifest, 'accepted_arithmetic_result': record(old / 'arithmetic-result.json'),
            'independent_terminal_review': record(old / 'independent-terminal-review.json'),
            'runtime_context': record(old / 'runtime-context.json'), 'build_review': record(old / 'build-review.json'),
            'owner': record(old / 'execute-gpu-stock/owner.json'), 'consumed_inputs': record(old / 'execute-gpu-stock/inputs-before.json'),
            'actual': record(old / 'execute-gpu-stock/actual.bin'), 'gpu': record(products / 'gpu'),
            'stock_module': record(products / 'stock.spv'), 'original_input': record(old / 'input.bin'),
            'original_expected': record(old / 'expected.bin'), 'all737_static_inputs_and17_terminal_owners_joined': True,
            'scope': 'Exact accepted complete stock arithmetic/source/build/producer/corpus/provider receipt reused. No fresh stock compile or dispatch is claimed.'}

def pin_context():
    verify_static()
    verify_host_reuse()
    verify_stock_reuse()
    verify_reference_acceptances()
    corpus_checks()
    assert not T.exists() and not (W / 'operator-context.json').exists()
    assert not subprocess.check_output(['git', 'status', '--porcelain=v1'], cwd=R)
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=R, text=True).strip()
    assert revision == P['runner_preparation_revision']
    actions = json.loads((W / 'actions.json').read_text())
    assert all(action['revision'] == revision for action in actions.values())
    write(W / 'operator-context.json', {'revision': revision, 'actions': record(W / 'actions.json'),
          'source_freeze': record(W / 'source-freeze.json'), 'preparation': record(W / 'preparation.json')})
    print('Pinned original multiply64 arithmetic preparation; no compiler/provider launched.', flush=True)

def dependencies(name):
    text = (T / name).read_text().replace('\\\n', '')
    return [Path(path).resolve(strict=True) for path in shlex.split(text.split(':', 1)[1])]

def module_checks():
    from interface_signature import interface_signature, function_body, normalized_function, multiply_shape
    texts = {name: (T / (name + '.spvasm')).read_text() for name in ('stock', 'candidate')}
    assert (T / 'stock.spv').read_bytes() == Path(P['original_stock_module']).read_bytes()
    unchanged = ('PB32Divide', 'PB32Unpack', 'PB32Round')
    assert all(normalized_function(texts['stock'], name) == normalized_function(texts['candidate'], name) for name in unchanged)
    results = {}
    for name, text in texts.items():
        caps = re.findall(r'OpCapability (\S+)', text)
        assert sorted(caps) == (['Shader'] if name == 'stock' else ['Int64', 'Shader'])
        assert 'OpTypeFloat' not in text
        assert not re.search(r'\b(?:RelaxedPrecision|FPFastMathMode|FPFastMathDefault|OpFmaKHR)\b', text)
        widths = dict(re.findall(r'(%\S+) = OpTypeInt (\d+) [01]', text))
        assert set(widths.values()) == ({'32'} if name == 'stock' else {'32', '64'})
        function = function_body(text, 'PB32Divide')
        divisions = [int(widths[kind]) for kind in re.findall(r'= OpUDiv (%\S+) ', function)]
        assert divisions == [32] and function.count('OpLoopMerge') == 1
        shape = None
        if name == 'candidate':
            assert re.findall(r'OpTypeInt 64 ([01])', text) == ['0']
            assert all(widths[kind] == '32' for kind in re.findall(r'= Op[US](?:Div|Mod|Rem) (%\S+) ', text))
            shape = multiply_shape(text)
        assert re.findall(r'OpExecutionMode %\S+ (.*)', text) == ['LocalSize 64 1 1']
        assert text.count('OpControlBarrier') == text.count('OpMemoryBarrier') == 0
        results[name] = {'module': record(T / (name + '.spv')), 'assembly': record(T / (name + '.spvasm')),
                         'divide_widths': divisions, 'divide_loops': function.count('OpLoopMerge'),
                         'capabilities': sorted(caps), 'integer_widths': sorted(set(widths.values())),
                         'multiply_graph': shape, 'original_divide_unpack_round_normalized_bodies_equal': True}
    assert interface_signature(texts['stock']) == interface_signature(texts['candidate'])
    return results

def build_review():
    seals = {label: verify_output_seal(label) for label in P['build_labels']}
    pins = json.loads((W / 'source-freeze.json').read_text())['pins']
    for depfile in ('host.d',):
        for path in dependencies(depfile):
            assert str(path) in pins and record(path) == pins[str(path)], path
    assert W / 'candidate/src/sirius/kernels/portable_binary32.h' in dependencies('host.d')
    assert W / 'stock/src/sirius/kernels/portable_binary32.h' not in dependencies('host.d')
    assert '-DSIRIUS_PB32_MULTIPLY64=1' in P['recipes']['compile-host'][0]
    assert '-DSIRIUS_PB32_MULTIPLY64=1' in P['recipes']['compile-shader-candidate'][0]
    needed = json.loads(Path(P['original_symbols']).read_text())['required_strong_sirius_definitions']
    strong = lambda text: {row.split()[2] for row in text.splitlines()
                           if len(row.split()) == 3 and row.split()[1] in ('T', 'D', 'B', 'R')}
    symbols = {}
    for name, binary in [('device-object', Path(P['device_object'])), ('gpu', T / 'gpu')]:
        text = subprocess.check_output(['/usr/bin/nm', '-g', '--defined-only', str(binary)], text=True)
        assert set(needed) <= strong(text)
        symbols[name] = text
    link_map = (T / 'gpu.map').read_text()
    assert 'libsirius_backend.a(vulkan_device.cpp.o)' not in link_map
    argv = P['reused_gpu_link_recipe']
    assert argv.index(P['device_object']) < argv.index(str(R / 'bin/linux-gcc/src/sirius/backend/libsirius_backend.a'))
    assert (T / 'resources/tests/backend/portable_binary32_reference.bin').read_bytes() == (W / 'portable_binary32_reference.bin').read_bytes()
    return {'passed': True, 'modules': module_checks(), 'output_seals': seals,
            'original_corpus': corpus_checks(),
            'same_admitted_host_reuse': verify_host_reuse(), 'required_strong_definitions': needed,
            'stock_arithmetic_receipt': verify_stock_reuse(), 'compiler_and_jit_receipts': verify_reference_acceptances(),
            'nm_outputs': symbols, 'private_object_before_original_archive': True,
            'original_vulkan_archive_member_extracted': False,
            'host_uses_candidate_header_and_explicit_macro': True,
            'probe_source_and_all_original_assertions_unchanged': True}

def phase(label):
    verify_static()
    context = json.loads((W / 'operator-context.json').read_text())
    assert context['actions'] == record(W / 'actions.json')
    assert context['source_freeze'] == record(W / 'source-freeze.json')
    assert context['preparation'] == record(W / 'preparation.json')
    assert label in P['outputs']
    for prior in P['dependencies'][label]:
        verify_output_seal(prior)
    assert all(not Path(path).exists() for path in P['outputs'][label])
    if label == 'prepare-build':
        assert not T.exists()
        (T / 'resources/tests/backend').mkdir(parents=True)
        (T / 'resources/tests/backend/portable_binary32_reference.bin').write_bytes((W / 'portable_binary32_reference.bin').read_bytes())
        for destination, original in P['reused_files'].items():
            shutil.copyfile(original, destination)
        (T / 'gpu').chmod(0o755)
    elif label == 'verify-build':
        write(W / 'build-review.json', build_review())
    else:
        for argv in P['recipes'][label]:
            subprocess.run(argv, cwd=R, check=True)
    write(W / label / 'output-seal.json', {'label': label, 'recipe': P['recipes'][label],
          'source_freeze': record(W / 'source-freeze.json'), 'operator_context': record(W / 'operator-context.json'),
          'inputs_before': record(W / label / 'inputs-before.json'),
          'outputs': {name: record(name) for name in P['outputs'][label]}})

def pin_runtime():
    verify_static()
    assert not (W / 'runtime-context.json').exists()
    verify_output_seal('verify-build')
    build = json.loads((W / 'build-review.json').read_text())
    assert build['passed'] and build['modules'] == module_checks()
    write(W / 'runtime-context.json', {'build_review': record(W / 'build-review.json'),
          'all_build_output_seals': {label: verify_output_seal(label) for label in P['outputs']},
          'gpu': record(T / 'gpu'), 'host': record(T / 'host'), 'modules': build['modules'],
          'original_input': record(W / 'input.bin'), 'original_expected': record(W / 'expected.bin'),
          'same_admitted_host_reuse': verify_host_reuse(),
          'scope': 'Original integer-only primitive corpus; no retained-program/table or normal-sum runtime acceptance.'})
    print('Original arithmetic runtime inputs pinned; no device launched.', flush=True)

def provider_mappings(owner):
    expected = set(P['provider_identities'])
    observed = set()
    births = {(row['pid'], row['saved_start_ticks']) for row in owner['observed_births']}
    for row in owner['provider_mappings']:
        assert (row['pid'], row['start_ticks']) in births
        assert row['owned_session'] == owner['child_pid']
        for line in row['maps']:
            fields = line.split(None, 5)
            assert len(fields) == 6 and fields[5] in expected
            identity = P['provider_identities'][fields[5]]
            device = [int(part, 16) for part in fields[3].split(':')]
            assert device == identity['device'] and int(fields[4]) == identity['inode']
            observed.add(fields[5])
    assert observed == expected

def join_runtime_inputs(label, owner):
    action = json.loads((W / 'actions.json').read_text())[label]
    folder = W / label
    assert owner['argv'] == [arg.replace('{folder}', str(folder)) for arg in action['argv']]
    before = json.loads((folder / 'inputs-before.json').read_text())
    for name in action['frozen_inputs']:
        path = Path(name)
        key = str(path.relative_to(R)) if path.is_relative_to(R) else str(path)
        assert before[key] == record(path), ('numerical consumed input differs from current sealed input', label, name)

def review():
    verify_static()
    context = json.loads((W / 'operator-context.json').read_text())
    assert context['actions'] == record(W / 'actions.json')
    assert context['preparation'] == record(W / 'preparation.json')
    runtime = json.loads((W / 'runtime-context.json').read_text())
    assert runtime['build_review'] == record(W / 'build-review.json')
    assert runtime['all_build_output_seals'] == {label: verify_output_seal(label) for label in P['outputs']}
    assert runtime['gpu'] == record(T / 'gpu') and runtime['host'] == record(T / 'host')
    assert runtime['modules'] == module_checks()
    assert runtime['original_input'] == record(W / 'input.bin') and runtime['original_expected'] == record(W / 'expected.bin')
    assert runtime['same_admitted_host_reuse'] == verify_host_reuse()
    for label in ('execute-host', 'execute-gpu-candidate'):
        owner = accepted(label)
        join_runtime_inputs(label, owner)
        if label.startswith('execute-gpu-'):
            provider_mappings(owner)
            assert (W / label / 'actual.bin').read_bytes() == (W / 'expected.bin').read_bytes()
            rows = [json.loads(line) for line in (W / label / 'stdout.log').read_text().splitlines() if line.startswith('{')]
            assert len(rows) == 1 and rows[0]['records'] == 174268
            assert rows[0]['complete_words'] == 697072 and rows[0]['mismatches'] == 0
            assert rows[0]['driver'] == 'llvmpipe'
    xml = ET.parse(W / 'execute-host/gtest.xml').getroot()
    assert int(xml.get('tests', '-1')) == 1 and not any(xml.findall('.//' + tag) for tag in ('failure', 'error', 'skipped'))
    cases = xml.findall('.//testcase')
    assert len(cases) == 1 and cases[0].get('classname') + '.' + cases[0].get('name') == 'PortableBinary32.MatchesIndependentExactCorpus'
    assert cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
    properties = {node.get('name'): node.get('value') for node in cases[0].findall('./properties/property')}
    assert properties['independent_exact_records_checked'] == '174268' and properties['exact_raw_word_mismatches'] == '0'
    write(W / 'arithmetic-result.json', {'schema': 'sirius-multiply64-original-arithmetic-result-v1',
          'passed': True, 'original_records': 174268, 'original_operation_count': 14,
          'candidate_macro': P['candidate_macro'], 'candidate_header': record(W / 'candidate/src/sirius/kernels/portable_binary32.h'),
          'candidate_shader': record(T / 'candidate.spv'), 'modules': module_checks(),
          'stock_result_reused': True, 'stock_arithmetic_receipt': verify_stock_reuse(),
          'compiler_and_jit_receipts': verify_reference_acceptances(), 'build_review': record(W / 'build-review.json'),
          'complete_output_words_per_provider': 697072, 'host_explicit64_and_both_full_gpu_outputs_equal_original_expected': True,
          'runtime_context': record(W / 'runtime-context.json'), 'original_corpus': corpus_checks(),
          'owners': {label: record(W / label / 'owner.json') for label in [*P['outputs'], 'execute-host', 'execute-gpu-candidate']},
          'scope': 'Full unchanged original primitive corpus; fresh candidate host/provider outputs and exact historical stock receipt. No retained programs, native-FP64 residual profile, integrated controls, performance, adoption or qualification verdict.'})
    print('Original174268 candidate host/provider outputs accepted; exact complete stock receipt reused.', flush=True)

if __name__ == '__main__':
    if sys.argv[1] == 'phase':
        phase(sys.argv[2])
    else:
        {'pin-context': pin_context, 'pin-runtime': pin_runtime, 'review': review}[sys.argv[1]]()
