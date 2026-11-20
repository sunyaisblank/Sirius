"""Run the original frozen C8 once, with its original reader and owner limits."""

from pathlib import Path
from types import SimpleNamespace
import hashlib
import json
import shutil
import subprocess
import sys
import zipfile

from prepare_provider import QUALIFICATION, R
from resume_raw import GUARD, check_space, pin, read_json, restore, write_json
from run_value_abi import expect_pin


HERE = Path(__file__).resolve().parent
ORIGINAL = QUALIFICATION / '363aa0c2/native-cpu-closed-all-original-c8/MANIFEST.json'
BINARY_ARCHIVE = QUALIFICATION / 'f878a18f/native-cpu-normalized-arithmetic-frame-timeout/actual-products.zip'
BASIS = 'out/native-cpu-closed-all-cohort'
CASE = 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults'
ORIGINAL_READER = R / 'out/native-cpu-arithmetic-sharing/readers/native_evidence.py'
ORIGINAL_LAYOUT = ORIGINAL_READER.with_name('expected-readback-layout.json')
LIMITS = {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}


def sealed(bundle, records, name):
    data = bundle.read(name)
    if {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} != records[name]:
        raise RuntimeError(f'Original C8 archive member differs: {name}')
    return data


def run(scratch):
    check_space()
    if subprocess.check_output(['git', 'status', '--porcelain']).strip():
        raise RuntimeError('Run C8 from the clean candidate checkout')
    linked = read_json(scratch / 'linked-provider.json')
    provider = R / linked['provider']['path']
    expect_pin(provider, linked['provider'])
    original = read_json(ORIGINAL)['bundles']['source-and-evidence']
    expect_pin(R / original['path'], original)
    with zipfile.ZipFile(R / original['path']) as archive:
        records = original['members']
        action = json.loads(sealed(archive, records, BASIS + '/actions.json'))['control-8']
        original_result = json.loads(sealed(archive, records, BASIS + '/control-8-acceptance.json'))['record']
        original_xml = scratch / 'original-c8.xml'
        original_xml.write_bytes(sealed(archive, records, BASIS + '/control-8/gtest.xml'))
    assert action['case'] == CASE and action['limits'] == LIMITS and action['expected_exit'] == 0
    assert pin(ORIGINAL_READER)['sha256'] == '1b44f0678dbfc7ee455b3d49e170be884b294882896fd1ac3706ff5ebded20ab'

    executable = scratch / 'audited_backend_tests'
    original_executable = action['executable']
    expected = action['expected_identities'][original_executable]
    if executable.exists():
        expect_pin(executable, expected)
    else:
        with zipfile.ZipFile(BINARY_ARCHIVE) as archive:
            with archive.open('files/' + str(Path(original_executable).relative_to(R))) as stream:
                restore(stream, executable, expected)
        executable.chmod(0o755)
    icd = scratch / 'lvp_icd.json'
    write_json(icd, {'file_format_version': '1.0.0',
                     'ICD': {'library_path': str(provider), 'api_version': '1.4.305'}})

    work = scratch / 'c8-owner'
    work.mkdir()
    observer = HERE / 'observe_c8.py'
    commands = work / 'observe.gdb'
    commands.write_text('set pagination off\nset confirm off\nset debuginfod enabled off\n'
                        'set auto-load off\nset startup-with-shell off\nset disable-randomization off\n'
                        f'source {observer}\nrun\n')
    shutil.copyfile(GUARD / 'run_owned.py', work / 'run_owned.py')
    (work / 'owned_guard.py').symlink_to(GUARD / 'owned_guard.py')
    folder = work / 'control-8'
    argv = ['/usr/bin/gdb', '-nx', '-nh', '-q', '-batch', '-x', str(commands), '--args',
            str(executable), '--gtest_filter=' + CASE, '--gtest_color=no']
    runtime_paths = [Path(line.split('=> ', 1)[1].split(' ', 1)[0]) for path in (provider, executable)
                     for line in subprocess.check_output(['/usr/bin/ldd', str(path)], text=True).splitlines()
                     if '=> /' in line]
    frozen = [executable, provider, icd, commands, observer, ORIGINAL_READER, ORIGINAL_LAYOUT,
              ORIGINAL, scratch / 'linked-provider.json', *runtime_paths]
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    write_json(work / 'actions.json', {'control-8': {
        'revision': revision, 'kind': 'test', 'expected_exit': 0, 'software_provider': True,
        'runtime_icd': str(icd), 'runtime_provider': str(provider), 'disable_software_shader_cache': True,
        'limits': LIMITS, 'argv': argv, 'frozen_inputs': sorted({str(path) for path in frozen}),
        'build_environment': {'SIMD16_C8_TASK': str(scratch), 'SIMD16_C8_OUTPUT': str(folder),
                              'PYTHONDONTWRITEBYTECODE': '1', 'TMPDIR': str(scratch)},
        'scope': 'Original candidate C8, unchanged executable, inputs and owner limits with passive SIMD16 observation',
    }})
    subprocess.run([sys.executable, '-B', str(work / 'run_owned.py'), 'control-8'], check=True)
    owner = read_json(folder / 'owner.json')
    assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['provider_mappings'] and not owner['remaining_owned_processes'] and not owner['cleanup_errors']

    sys.path.insert(0, str(ORIGINAL_READER.parent))
    import native_evidence as reader
    modules = original_result['selected_modules']
    names = {stage: value['array'] for stage, value in modules.items()}
    reader.selected_names = lambda context, arm, wide=False: names
    original_props = reader.xml_properties(original_xml, 'control-8')
    features = {key: original_props[key] for key in ('fma_fp32_enabled', 'preserves_fp32_signed_zero_inf_nan')}
    reference_context = SimpleNamespace(value={'native_identity': {
        **original_result['identity'], **features}},
        arms={'candidate': {'arrays': {value['array']: value for value in modules.values()}}})
    props = reader.xml_properties(folder / 'gtest.xml', 'control-8')
    result = reader.timestamp_result(props, 'candidate', reference_context)
    assert result['readbacks'] == original_result['readbacks'] and len(result['readbacks']) == 32

    observation = read_json(folder / 'native-observation.json')
    assert not observation['failures'] and observation['exit_code'] == 0
    packets = observation['packets']
    for key in ('transport:3:1', 'transport:12:1', 'endpoint:3:1', 'endpoint:12:8',
                'dopri_phase:24:16', 'dopri_phase:24:8'):
        assert packets[key]['handled'], key
    wide = packets['dopri_phase:24:16']
    remainder = packets['dopri_phase:24:8']
    assert wide['x'] == 0 and remainder['x'] == 16 and wide['width'] == remainder['width'] == 16
    assert wide['callback_restored_mxcsr'] == wide['callback_saved_mxcsr']
    assert observation['wide_dp']['thread'] == wide['thread'] and 'return_mxcsr' in observation['wide_dp']
    write_json(scratch / 'c8-acceptance.json', {'source': revision, 'provider': linked['provider'],
        'original_test': expected, 'owner': pin(folder / 'owner.json'), 'xml': pin(folder / 'gtest.xml'),
        'observation': pin(folder / 'native-observation.json'), 'record': result, 'native': observation})


if __name__ == '__main__':
    location, = (Path(arg).resolve() for arg in sys.argv[1:])
    if location == R / 'out' or not location.is_relative_to(R / 'out'):
        raise RuntimeError('Use the working directory for this provider build under out/')
    run(location)
