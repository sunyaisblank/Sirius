"""Finish the original finite controls after the sole ABBA and Decimal pass."""

from pathlib import Path
from types import SimpleNamespace
import json
import shutil
import subprocess
import sys
import zipfile

from run_original_c8 import (GUARD, LIMITS, ORIGINAL_LAYOUT, ORIGINAL_READER, QUALIFICATION,
                             R, check_space, expect_pin, pin, read_json, sealed, write_json)


AX_FRAME = QUALIFICATION / '0fa8215e/native-cpu-closed-all-original-frame-attempt/MANIFEST.json'
AX_ROOT = 'out/native-cpu-closed-all-cohort'


def run(scratch):
    check_space()
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    if subprocess.check_output(['git', 'status', '--porcelain']).strip():
        raise RuntimeError('Run the original controls from the clean C8 candidate checkout')
    c8 = read_json(scratch / 'c8-acceptance.json')
    abba = read_json(scratch / 'abba/benefit.json')
    decimal = read_json(scratch / 'abba/decimal-agreement.json')
    assert c8['source'] == abba['source'] == revision
    assert abba['all22_numeric_pass'] and decimal['agrees_on_all22'] and decimal['all22_pass']
    original = read_json(AX_FRAME)['bundle']
    expect_pin(R / original['path'], original)
    with zipfile.ZipFile(R / original['path']) as archive:
        inputs = original['members']
        original_actions = json.loads(sealed(archive, inputs, AX_ROOT + '/controls/runtime/actions.json'))
        original_result = json.loads(sealed(archive, inputs, AX_ROOT + '/controls/readers/controls-1-8-acceptance.json'))
    original_routes = original_result['remaining7']['control-1']['complete_native_routes']['routes']

    provider = R / c8['provider']['path']
    executable = scratch / 'audited_backend_tests'
    icd = scratch / 'lvp_icd.json'
    expect_pin(provider, c8['provider'])
    expect_pin(executable, c8['original_test'])
    work = scratch / 'remaining-controls'
    work.mkdir()
    shutil.copyfile(GUARD / 'run_owned.py', work / 'run_owned.py')
    (work / 'owned_guard.py').symlink_to(GUARD / 'owned_guard.py')
    sys.path.insert(0, str(ORIGINAL_READER.parent))
    import native_evidence as reader

    frozen = [provider, executable, icd, ORIGINAL_READER, ORIGINAL_LAYOUT,
              Path(__file__).resolve(), AX_FRAME, scratch / 'c8-acceptance.json',
              scratch / 'abba/benefit.json', scratch / 'abba/decimal-agreement.json']
    actions = {}
    for i in range(1, 8):
        label = 'control-' + str(i)
        old = original_actions[label]
        assert old['case'] == reader.CASES[label] and old['limits'] == LIMITS and old['expected_exit'] == 0
        actions[label] = {
            'revision': revision, 'kind': 'test', 'expected_exit': 0, 'software_provider': True,
            'runtime_icd': str(icd), 'runtime_provider': str(provider), 'disable_software_shader_cache': True,
            'limits': LIMITS, 'frozen_inputs': list(map(str, frozen)),
            'argv': [str(executable), '--gtest_filter=' + old['case'], '--gtest_color=no'],
            'scope': 'Original finite control, unchanged executable, inputs, cache and resource limits',
            'build_environment': {'TMPDIR': str(scratch), 'PYTHONDONTWRITEBYTECODE': '1'},
        }
    write_json(work / 'actions.json', actions)
    results = {}
    for label, action in actions.items():
        subprocess.run([sys.executable, '-B', str(work / 'run_owned.py'), label], check=True)
        folder = work / label
        owner = read_json(folder / 'owner.json')
        assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert owner['provider_mappings'] and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
        props = reader.xml_properties(folder / 'gtest.xml', label)
        result = {'case': original_actions[label]['case'], 'wall_seconds': owner['wall_seconds'],
                  'peak_sampled_rss_kib': owner['peak_sampled_rss_kib'],
                  'owner': pin(folder / 'owner.json'), 'xml': pin(folder / 'gtest.xml')}
        if label == 'control-1':
            routes = {item['wide']: item['selected_arrays'] for item in original_routes}
            arrays = {name: {'sha256': item['selected_stage_sha256'][stage]}
                      for item in original_routes for stage, name in item['selected_arrays'].items()}
            c8_props = reader.xml_properties(scratch / 'c8-owner/control-8/gtest.xml', 'control-8')
            context = SimpleNamespace(value={'native_identity': {'fma_fp32_enabled': c8_props['fma_fp32_enabled']}},
                                      arms={'candidate': {'arrays': arrays}})
            reader.selected_names = lambda ctx, arm, wide=False: routes[wide]
            facts = reader.factored_routes({'props': props, 'folder': folder, 'arm': 'candidate'}, context)
            assert facts['routes_count'] == 2 and facts['complete_readbacks'] == 28
            result['factored_routes'] = facts
        results[label] = result
    write_json(scratch / 'all8-acceptance.json', {'source': revision, 'provider': c8['provider'],
        'all8_original_controls_pass': True, 'original_c8': pin(scratch / 'c8-acceptance.json'),
        'original_abba': pin(scratch / 'abba/benefit.json'), 'decimal': pin(scratch / 'abba/decimal-agreement.json'),
        'remaining7': results})


if __name__ == '__main__':
    location, = (Path(arg).resolve() for arg in sys.argv[1:])
    if location == R / 'out' or not location.is_relative_to(R / 'out'):
        raise RuntimeError('Use the working directory for this provider build under out/')
    run(location)
