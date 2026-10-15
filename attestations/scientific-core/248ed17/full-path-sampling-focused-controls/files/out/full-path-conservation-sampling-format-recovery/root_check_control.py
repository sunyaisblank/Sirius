"""Root passive terminal/XML/input check; never starts or retries a control."""
from pathlib import Path
import hashlib
import json
import os
import sys
import xml.etree.ElementTree as ET

R = Path.cwd()
Q = Path(__file__).resolve().parent
C = Q / 'controls'
case, chunk = sys.argv[1:]
a = json.loads((C / 'actions.json').read_text())[case]
d = C / case
o = json.loads((d / 'owner.json').read_text())


def pin(path):
    path = Path(path).resolve(strict=True)
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return {'bytes': path.stat().st_size, 'sha256': h.hexdigest(), 'resolved_path': str(path)}


def ticks(pid):
    try:
        return Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
    except (FileNotFoundError, ProcessLookupError):
        return None


assert o['source'] == {'revision': a['revision'], 'status': ''}
assert o['case'] == case and o['argv'] == a['argv'] and o['limits'] == a['limits']
assert o['status'] == 'terminal' and o['passed'] is True and o['child_exit'] == 0
assert o['stop_reason'] is None and o['cleanup_errors'] == o['remaining_owned_processes'] == []
assert o['source_unchanged'] is o['whole_inputs_unchanged'] is True
assert o['wall_seconds'] < a['limits']['seconds']
assert o['peak_sampled_rss_kib'] <= a['limits']['rss_kib']
assert o['minimum_sampled_available_kib'] >= a['limits']['reserve_kib']
assert o['boot_id'] == Path('/proc/sys/kernel/random/boot_id').read_text().strip()
before = json.loads((d / 'inputs-before.json').read_text())
assert before == json.loads((d / 'inputs-after.json').read_text())
for name, wanted in before.items():
    assert pin(R / name) == wanted, name
births = []
for row in o['observed_births']:
    current = ticks(row['pid'])
    assert row['absent'] is True and current != str(row['saved_start_ticks'])
    births.append({'pid': row['pid'], 'saved_start_ticks': str(row['saved_start_ticks']),
                   'current_start_ticks': current, 'absent': True})
assert births and all(row['wait_status'] == 0 for row in o['adopted_children_reaped'])
x = ET.parse(d / 'gtest.xml').getroot()
suites, tests = x.findall('testsuite'), x.findall('.//testcase')
assert len(suites) == len(tests) == 1
w = a['root_required_xml']
for parent in [x, suites[0]]:
    assert int(parent.attrib['tests']) == 1
    assert all(int(parent.attrib.get(name, '0')) == 0 for name in ['failures', 'errors', 'disabled', 'skipped'])
assert suites[0].attrib['name'] == w['suite']
assert all(tests[0].attrib[name] == w[name] for name in ['name', 'status', 'result'])
assert tests[0].attrib['classname'] == w['suite']
assert not any(x.findall('.//' + tag) for tag in ['failure', 'error', 'skipped'])
env = {'GTEST_OUTPUT': 'xml:' + str(d / 'gtest.xml')}
if a['software_provider']:
    env.update(VK_DRIVER_FILES='/usr/share/vulkan/icd.d/lvp_icd.json', SIRIUS_VULKAN_DEVICE='0')
    assert a['disable_software_shader_cache'] is False and o['provider_mappings']
    provider = Path('/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so')
    st = provider.stat()
    for mapping in o['provider_mappings']:
        assert mapping['owned_session'] == o['child_pid']
        assert any(row['pid'] == mapping['pid'] and row['saved_start_ticks'] == str(mapping['start_ticks']) for row in births)
        assert mapping['maps']
        for line in mapping['maps']:
            fields = line.split()
            assert fields[5] == str(provider) and int(fields[4]) == st.st_ino
            assert tuple(int(v, 16) for v in fields[3].split(':')) == (os.major(st.st_dev), os.minor(st.st_dev))
else:
    assert o['provider_mappings'] == []
assert o['declared_output_environment'] == env
result = {'schema': 'sirius-full-path-sampling-root-control-terminal-check-v1',
          'passed': True, 'source_revision': a['revision'], 'case': case,
          'actual_root_command': {'executed_once': True, 'closed': True, 'exit_code': 0, 'terminal_tool_chunk': chunk},
          'current_whole_input_rows': len(before), 'before_after_rows_checked': 2 * len(before),
          'fresh_observed_births_absent': births, 'wall_seconds': o['wall_seconds'],
          'peak_sampled_rss_kib': o['peak_sampled_rss_kib'], 'exact_one_case_XML_zero_FESD': True,
          'actual_LVP_mapping_required_and_bound': a['software_provider'],
          'default_shader_cache': True, 'properties': {p.attrib['name']: p.attrib['value'] for p in tests[0].findall('./properties/property')},
          'evidence': {name: pin(d / name) for name in ['owner.json', 'inputs-before.json', 'inputs-after.json', 'gtest.xml', 'stdout.log', 'stderr.log', 'samples.jsonl']},
          'scope': 'One unchanged original focused control only; prospective outer wall is diagnostic. No full scientific/domain/native/release qualification.'}
output = Q / ('root-control-' + case + '.json')
with output.open('x') as stream:
    json.dump(result, stream, indent=2, sort_keys=True)
    stream.write('\n')
print(json.dumps({'passed': True, 'case': case, 'receipt': pin(output), 'input_rows': len(before), 'births': len(births), 'wall_seconds': o['wall_seconds']}))
