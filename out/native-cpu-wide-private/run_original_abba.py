"""Run the sole original A1/B1/B2/A2 sequence after current C8 acceptance."""

from pathlib import Path
from types import SimpleNamespace
import json
import shutil
import subprocess
import sys
import zipfile

from run_original_c8 import (CASE, GUARD, LIMITS, ORIGINAL_LAYOUT, ORIGINAL_READER, QUALIFICATION,
                             R, check_space, expect_pin, pin, read_json, restore, sealed, write_json)


AX_BENEFIT = QUALIFICATION / 'f9428502/native-cpu-closed-all-original-benefit/MANIFEST.json'
AX_DIR = 'out/native-cpu-closed-all-cohort/benefit'
BASELINE_BUNDLE = QUALIFICATION / '08add7e5/native-original-input-ownership-finite-outcome/actual-products-03.zip'
DECIMAL = R / 'out/native-original-input-ownership/readers/decimal_check.py'
NUMERIC = ORIGINAL_READER.with_name('numeric_gates.py')
ORDER = ('a1', 'b1', 'b2', 'a2')


def run(scratch):
    check_space()
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    if subprocess.check_output(['git', 'status', '--porcelain']).strip():
        raise RuntimeError('Run ABBA from the clean C8 candidate checkout')
    c8 = read_json(scratch / 'c8-acceptance.json')
    assert c8['source'] == revision and len(c8['record']['readbacks']) == 32
    c8_owner_path = scratch / 'c8-owner/control-8/owner.json'
    expect_pin(c8_owner_path, c8['owner'])
    c8_owner = read_json(c8_owner_path)
    assert c8_owner['passed']
    original = read_json(AX_BENEFIT)['bundle']
    expect_pin(R / original['path'], original)
    with zipfile.ZipFile(R / original['path']) as archive:
        records = original['members']
        original_actions = json.loads(sealed(archive, records, AX_DIR + '/runtime/actions.json'))
        reference = json.loads(sealed(archive, records, AX_DIR + '/readers/benefit.json'))['records']
        xml = {'baseline': 'a1', 'candidate': 'b1'}
        for arm, slot in xml.items():
            (scratch / f'original-abba-{arm}.xml').write_bytes(sealed(archive, records, f'{AX_DIR}/runtime/{slot}/gtest.xml'))

    old_a = original_actions['a1']
    baseline = scratch / 'baseline-libvulkan_lvp.so'
    baseline_pin = old_a['expected_identities'][old_a['runtime_provider']]
    with zipfile.ZipFile(BASELINE_BUNDLE) as archive, archive.open('sha256/' + baseline_pin['sha256']) as stream:
        restore(stream, baseline, baseline_pin)
    baseline_icd = scratch / 'baseline-lvp_icd.json'
    write_json(baseline_icd, {'file_format_version': '1.0.0',
                             'ICD': {'library_path': str(baseline), 'api_version': '1.4.305'}})
    candidate = R / c8['provider']['path']
    expect_pin(candidate, c8['provider'])
    executable = scratch / 'audited_backend_tests'
    expect_pin(executable, c8['original_test'])
    provider = {'baseline': (baseline, baseline_icd), 'candidate': (candidate, scratch / 'lvp_icd.json')}
    runtime = scratch / 'abba'
    runtime.mkdir()
    shutil.copyfile(GUARD / 'run_owned.py', runtime / 'run_owned.py')
    (runtime / 'owned_guard.py').symlink_to(GUARD / 'owned_guard.py')
    sys.path.insert(0, str(ORIGINAL_READER.parent))
    import native_evidence as reader
    from numeric_gates import evaluate_validated_records

    contexts, names = {}, {}
    for arm, slot in xml.items():
        props = reader.xml_properties(scratch / f'original-abba-{arm}.xml', slot)
        saved = reference[slot]
        modules = saved['selected_modules']
        names[arm] = {stage: row['array'] for stage, row in modules.items()}
        hardware = {**saved['identity'],
                    **{key: props[key] for key in ('fma_fp32_enabled', 'preserves_fp32_signed_zero_inf_nan')}}
        contexts[arm] = SimpleNamespace(value={'native_identity': hardware},
            arms={arm: {'arrays': {row['array']: row for row in modules.values()}}})
    reader.selected_names = lambda context, arm, wide=False: names[arm]

    assert pin(NUMERIC)['sha256'] == '051b1d0a83c5669099ca00b0496302d662d001b0f341d291a7cf89f68e3abdf1'
    assert pin(DECIMAL)['sha256'] == 'd27a3e94775edbdba42b48569d0e09a794b720fc13bc94d92787f9c62b1f305b'
    shared = [executable, ORIGINAL_READER, ORIGINAL_LAYOUT, AX_BENEFIT, c8_owner_path,
              scratch / 'c8-acceptance.json', DECIMAL, NUMERIC, Path('/usr/lib/x86_64-linux-gnu/libvulkan.so.1')]
    actions = {}
    for slot in ORDER:
        old = original_actions[slot]
        assert old['case'] == CASE and old['limits'] == LIMITS and old['expected_exit'] == 0
        arm = old['arm']
        library, icd = provider[arm]
        actions[slot] = {'revision': revision, 'kind': 'test', 'expected_exit': 0, 'software_provider': True,
            'disable_software_shader_cache': True, 'runtime_provider': str(library), 'runtime_icd': str(icd),
            'limits': LIMITS, 'frozen_inputs': [str(p) for p in (*shared, library, icd)],
            'argv': [str(executable), '--gtest_filter=' + CASE, '--gtest_color=no'], 'arm': arm,
            'scope': 'Original sole ABBA timing slot, unchanged executable, inputs, cache and resource limits',
            'build_environment': {'TMPDIR': str(scratch), 'PYTHONDONTWRITEBYTECODE': '1'}}
    write_json(runtime / 'actions.json', actions)
    accepted, owners = {}, {}
    last = c8_owner['child_start_ticks']
    for slot in ORDER:
        subprocess.run([sys.executable, '-B', str(runtime / 'run_owned.py'), slot], check=True)
        folder = runtime / slot
        owner = read_json(folder / 'owner.json')
        assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert owner['provider_mappings'] and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
        assert owner['boot_id'] == c8_owner['boot_id'] and owner['child_start_ticks'] > last
        last = owner['child_start_ticks']
        owners[slot] = {'owner': pin(folder / 'owner.json'), 'xml': pin(folder / 'gtest.xml'),
                        'wall_seconds': owner['wall_seconds'], 'peak_sampled_rss_kib': owner['peak_sampled_rss_kib']}
        arm = actions[slot]['arm']
        accepted[slot] = reader.timestamp_result(reader.xml_properties(folder / 'gtest.xml', slot), arm, contexts[arm])
        assert accepted[slot]['readbacks'] == c8['record']['readbacks']
    numeric = evaluate_validated_records(accepted, {'observed_instance_api_minimum_improvement_percent': 0})
    summary = {'source': revision, 'fixed_order': ORDER, 'audit_passed': True,
        'benefit_passed': numeric['all22_numeric_pass'], 'owners': owners, **numeric}
    write_json(runtime / 'benefit.json', summary)
    source = DECIMAL.read_text()
    anchor = "R=Path.cwd();W=Path(__file__).resolve().parent;RUN=R/'out/native-original-input-ownership/runtime'"
    assert source.count(anchor) == 1
    bound = source.replace(anchor, f'W=RUN=Path({str(runtime)!r})')
    exec(compile(bound, str(DECIMAL), 'exec'), {})
    return summary


if __name__ == '__main__':
    location, = (Path(arg).resolve() for arg in sys.argv[1:])
    if location == R / 'out' or not location.is_relative_to(R / 'out'):
        raise RuntimeError('Use the working directory for this provider build under out/')
    sys.exit(0 if run(location)['benefit_passed'] else 1)
