"""Read existing evidence without rebuilding, rendering, or changing its files.

Run from the Sirius checkout. Output is JSON, not a qualification receipt.
"""
from pathlib import Path
import ast
import collections
from decimal import Decimal, localcontext
import hashlib
import json
import re
import struct
import subprocess
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path.cwd()
NATIVE = ROOT / 'out/native-cpu-slang-bridge'
RESULT = {}


def identity(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def check(data, expected):
    assert identity(data) == {key: expected[key] for key in ('bytes', 'sha256')}


def doc(path):
    return json.loads(Path(path).read_text())


def archive_bundle(folder, groups):
    manifest = doc(folder / 'MANIFEST.json')
    archives = {}
    for name, expected in manifest['archives'].items():
        check((folder / name).read_bytes(), expected)
        archives[name] = zipfile.ZipFile(folder / name)
    count = 0
    for group in groups:
        for name, expected in manifest[group].items():
            locations = ([expected['archive']] if 'archive' in expected else
                         [a for a, z in archives.items() if name in z.namelist()])
            if not locations:
                reused = manifest.get('reused_product_contents', {}).get(name)
                assert reused
                data = subprocess.check_output(['git', 'show',
                    reused['git_revision'] + ':' + reused['git_path']])
            else:
                assert len(locations) == 1, name
                data = archives[locations[0]].read(name)
            check(data, expected)
            count += 1
    return manifest, archives, count


folder = ROOT / 'attestations/qualification/536ac050/native-cpu-bridge-source-and-boundaries'
manifest, archives, count = archive_bundle(folder, ['members'])
RESULT['native_archive_integrity'] = {'archives': len(archives), 'members': count}
compiled = doc(NATIVE / 'batch-compile-result.json')
reused = doc(NATIVE / 'reused-row-products.json')
for name, expected in reused['restored_actual_products'].items():
    check((ROOT / name).read_bytes(), expected)
RESULT['reused_current_row_host_products_verified'] = len(reused['restored_actual_products'])
required = {'-std=c++2c', '-Werror', '-march=x86-64', '-fexcess-precision=standard',
            '-fno-fast-math', '-fno-associative-math', '-frounding-math', '-ffp-contract=off'}
for item in compiled['objects']:
    assert item['exit'] == 0 and required <= set(item['argv'])
    check((ROOT / item['object']['path']).read_bytes(), item['object'])
for name, expected in compiled['source_inputs'].items():
    check((ROOT / name).read_bytes(), expected)
build = doc(NATIVE / 'provider/provider-build-result.json')
library = ROOT / build['library']
check(library.read_bytes(), build['library_identity'])

# Reuse the inspected ELF parser definition, not the writer/runtime commands in
# the historical reader. Independently walk all actual object relocations below.
tree = ast.parse((NATIVE / 'read_linked_native.py').read_text())
definition = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == 'Elf')
namespace = {'struct': struct}
exec(compile(ast.Module(body=[definition], type_ignores=[]), '<inspected ELF parser>', 'exec'), namespace)
Elf = namespace['Elf']
linked = Elf(library)
functions = relocations_checked = 0
for item in build['native_objects']:
    path = ROOT / item['path']
    check(path.read_bytes(), item)
    obj = Elf(path)
    relocations = collections.defaultdict(list)
    for index, section in enumerate(obj.sections):
        if section[1] == 4:
            for offset in range(0, section[5], 24):
                location, info, _ = struct.unpack_from('<QQq', obj.section_data(index), offset)
                relocations[section[7]].append((location, info & 0xffffffff))
    for symbol in obj.symbols:
        if symbol['type'] != 2 or not symbol['size'] or not symbol['section']:
            continue
        if symbol['name'].startswith('_GLOBAL__sub_I_'):
            assert obj.function_data(symbol) == b'\xc3'
            assert not [s for s in linked.symbols if s['name'] == symbol['name']]
            continue
        final = linked.named(symbol['name'])
        assert final['size'] == symbol['size']
        left, right = bytearray(obj.function_data(symbol)), bytearray(linked.function_data(final))
        for location, kind in relocations[symbol['section']]:
            offset = location - symbol['value']
            if 0 <= offset < symbol['size']:
                assert kind in {2, 4, 42} and offset + 4 <= symbol['size']
                left[offset:offset+4] = right[offset:offset+4] = b'\0' * 4
                relocations_checked += 1
        assert left == right, symbol['name']
        functions += 1
assert functions == 64 and relocations_checked == 241
batch = doc(NATIVE / 'batch-generation.json')
for stage in ('Transport', 'Endpoint', 'DopriPhase'):
    item = batch['modules'][stage]['SPIRV']
    data = (ROOT / item['module']).read_bytes()
    check(data, item)
    assert linked.data.count(data) == 1
callback = subprocess.check_output(['objdump', '-d', '--disassemble=cs_exec_fn', str(library)], text=True)
for slug in ('transport', 'endpoint', 'dopri_phase'):
    assert '<sirius_native_' + slug + '_try_packet>' in callback
assert '<util_fpstate_get>' in callback and callback.count('<util_fpstate_set>') == 2
RESULT['native_current_products'] = {'objects': 3, 'linked_functions': functions,
    'declared_relocation_operands': relocations_checked, 'complete_module_spans': 3,
    'strict_compile_commands_and_sealed_inputs_verified': True,
    'successful_native_execution_proved': False}

guards = {}
for slug in ('transport', 'endpoint', 'dopri_phase'):
    p = NATIVE / 'ordinary' / slug
    owner = doc(p / 'owner.json')
    assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert all(x['absent'] for x in owner['observed_births'])
    before = doc(p / 'inputs-before.json')
    assert before == doc(p / 'inputs-after.json')
    for name, expected in before.items():
        check((ROOT / name).read_bytes(), expected)
    call, = doc(p / 'native-calls.json')
    boundary, = doc(p / 'compute-boundaries.json')
    verification = doc(p / 'capture-verification.json')
    assert call['returned'] and call['handled'] is False
    assert call['input_bytes'] == call['output_bytes'] == 4
    assert call['output_before'] == call['output_after'] == 'baaadbda'
    assert call['mxcsr'] & 0xe040 == 0
    assert boundary['returned'] and boundary['worker'] == call['thread']
    assert boundary['mxcsr_before'] == boundary['mxcsr_after']
    assert verification['passed'] and verification['callback_restoration_observed']
    assert verification['actual_entry']['worker_thread'] == call['thread']
    assert verification['actual_provider_sha256'] == build['library_identity']['sha256']
    guards[slug] = {'finite_negative_boundary_passed': True,
        'wall_seconds': owner['wall_seconds'], 'successful_arithmetic_exercised': False}
RESULT['negative_guards'] = guards

folder = ROOT / 'attestations/qualification/34ef39fb/native-cpu-arithmetic-sharing-finite-outcome'
manifest, archives, count = archive_bundle(folder, ['diagnostics', 'actual_products'])
z = archives['diagnostics-and-trial.zip']
prefix = 'out/native-cpu-arithmetic-sharing/'
owner = json.loads(z.read(prefix + 'cold-frame/frame/owner.json'))
assert not owner['passed'] and owner['status'] == 'terminal'
assert owner['stop_reason'] == 'diagnostic_timeout' and owner['wall_seconds'] >= 300
assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
assert z.read(prefix + 'cold-frame/frame/inputs-before.json') == z.read(prefix + 'cold-frame/frame/inputs-after.json')
log = z.read(prefix + 'cold-frame/frame/stdout.log').decode()
assert '2568' in log and '16384' in log
assert prefix + 'cold-frame/frame/gtest.xml' not in z.namelist()
xmls = [n for n in z.namelist() if n.endswith('/gtest.xml')]
for name in xmls:
    root = ET.fromstring(z.read(name))
    assert not any(int(node.attrib.get(k, 0)) for node in root.iter()
                   for k in ('failures', 'errors', 'disabled', 'skipped'))
RESULT['previous_candidate'] = {'members_verified': count, 'completed_zero_failure_xml_files': len(xmls),
    'frame_wall_seconds': owner['wall_seconds'], 'frame_passed': False,
    'last_logged_progress': [2568, 16384], 'completed_frame_xml_present': False,
    'scope': 'Historical candidate only; does not qualify the current native bridge.'}
benefit = json.loads(z.read(prefix + 'readers/benefit.json'))
decimal = json.loads(z.read(prefix + 'readers/decimal-agreement.json'))
assert benefit['fixed_order'] == ['a1', 'b1', 'b2', 'a2']
assert benefit['internal_repeats_per_slot'] == 3 and benefit['extra_retries_or_warmups'] == 0
assert len(decimal['checks']) == 22
with localcontext() as context:
    context.prec = 60
    for item in decimal['checks']:
        a, b = Decimal(item['baseline_ms']), Decimal(item['candidate_ms'])
        gain = 100 * (a-b) / a
        assert a > 0 and b > 0 and gain >= Decimal(item['minimum_percent'])
        assert abs(gain - Decimal(item['improvement_percent'])) < Decimal('1e-48')
reference = None
for label in ['control-8', 'a1', 'b1', 'b2', 'a2']:
    root = ET.fromstring(z.read(prefix + 'runtime/' + label + '/gtest.xml'))
    properties = {node.attrib['name']: node.attrib['value'] for node in root.iter('property')}
    readbacks = {k: dict(part.split('=', 1) for part in v.split(';'))
                 for k, v in properties.items() if '_readback_' in k}
    assert len(readbacks) == 32
    if reference is None:
        reference = readbacks
    assert readbacks == reference
    if label != 'control-8':
        assert readbacks == benefit['records'][label]['readbacks']
controls = json.loads(z.read(prefix + 'readers/controls-1-8-acceptance.json'))
assert controls['all8_original_controls_pass']
RESULT['previous_candidate'].update({'historical_sole_ABBA_Decimal_checks_recomputed': 22,
    'C8_and_all_four_slots_have_32_equal_complete_readback_digests': True,
    'all8_historical_control_XMLs_pass': True,
    'limit': 'Digest equality and archived original admission are reused; buffers are not 32 archived raw buffer files. No new workload or claim of fresh scientific qualification.'})

folder = ROOT / 'attestations/qualification/0109ab46/native-build-domains'
manifest = doc(folder / 'MANIFEST.json')
for name, expected in manifest['bundle_files'].items():
    check((folder / name).read_bytes(), expected)
native_gates = {}
for platform in ('macos', 'windows'):
    z = zipfile.ZipFile(folder / (platform + '-native-attestation.zip'))
    gates = [n for n in z.namelist() if n.endswith('native_build_gate.json')]
    assert len(gates) == 1
    gate = json.loads(z.read(gates[0]))
    assert gate['source']['revision'] == '0109ab462edade524667d20303507aa09adcbaf6'
    assert gate['source']['clean']
    ctest = gate['ctest']
    assert ctest['registered'] == 1152 and ctest['executed'] == ctest['selected'] == 9
    assert not any(ctest[k] for k in ('failures', 'errors', 'skipped'))
    native_gates[platform] = {'registered': 1152, 'executed': 9, 'full_estate_executed': False}
receipt = doc(folder / 'alignment-receipt.json')
assert receipt['admitted_domains'] == ['macos-native-build', 'windows-native-build']
RESULT['historical_public_build_domains'] = {'revision': receipt['source_revision'],
    'admitted': receipt['admitted_domains'], 'gates': native_gates,
    'admissible_for_current_private_revision': False}

pairs = {}
for name in ('kerr-png-reconstruction-pair', 'kerr-exr-reconstruction-pair', 'alcubierre-png-reconstruction-pair'):
    folder = ROOT / 'attestations/forensic-images/b739f92' / name
    result = doc(folder / 'result.json')
    assert result['evidence_kind'] == 'NEW_RECONSTRUCTION'
    assert not result['historical_tape_identity_claimed']
    assert result['comparison']['equal'] and result['comparison']['mismatched_samples'] == 0
    manifest = doc(folder / 'manifest.json')
    for filename, expected in manifest['files'].items():
        path = folder / filename
        if expected.get('type') == 'file':
            check(path.read_bytes(), expected)
    for item in result['inputs']:
        relative = str(Path(item['path']).relative_to(ROOT))
        check((folder / 'files' / relative).read_bytes(), item)
    pairs[name] = {'manifest_files_verified': len(manifest['files']),
        'compared_samples': result['comparison']['compared_samples'],
        'exact_mismatches': 0, 'original_July_provenance_established': False}
RESULT['reconstruction_pairs'] = pairs
controls = doc(ROOT / 'attestations/forensic-images/b739f92/kerr-png-reconstruction-pair/files/out/forensic-reconstruction-b739f92/comparator-controls/receipt.json')
assert controls['passed'] and not controls['historical_tape_identity_claim']
for command in controls['commands']:
    assert command['exit'] == command['expected_exit']
perturbed = [command for command in controls['commands']
             if any('perturbed.' in arg for arg in command['argv'])]
assert len(perturbed) == 3 and all(command['exit'] == 1 for command in perturbed)
RESULT['historical_comparator_controls'] = {'commands': len(controls['commands']),
    'perturbed_pixel_expected_failures': len(perturbed), 'recorded_exits_match': True,
    'fresh_execution_claimed': False}
RESULT['qualification_claim'] = False
print(json.dumps(RESULT, indent=2))
