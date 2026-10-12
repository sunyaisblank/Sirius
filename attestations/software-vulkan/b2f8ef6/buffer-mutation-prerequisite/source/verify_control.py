#!/usr/bin/env python3
"""Check the finite buffer-mutation prerequisite without admitting qualification."""
import hashlib
import json
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
BIN = ROOT / 'bin/linux-gcc/vulkan-buffer-mutation'
EXECUTED = 'b2f8ef65265caa11cc96b23dc6834dc1cd88d351'

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

reports = {}
for action in ('compile-vulkan', 'compile-test', 'link', 'control'):
    folder = WORK / action
    report = json.loads((folder / 'owner.json').read_text())
    assert report['source'] == {'revision': EXECUTED, 'status': ''}
    assert report['passed'] and report['child_exit'] == 0
    assert report['stop_reason'] is None and report['cleanup_errors'] == []
    assert report['remaining_owned_processes'] == []
    assert report['source_unchanged'] and report['whole_inputs_unchanged']
    assert all(item['absent'] for item in report['observed_births'])
    assert json.loads((folder / 'inputs-before.json').read_text()) == json.loads((folder / 'inputs-after.json').read_text())
    reports[action] = {k: report[k] for k in ('wall_seconds', 'peak_sampled_rss_kib', 'observed_births')}
report = json.loads((WORK / 'control/owner.json').read_text())
assert report['provider_mappings']
assert report['declared_output_environment']['VK_DRIVER_FILES'] == '/usr/share/vulkan/icd.d/lvp_icd.json'
assert report['declared_output_environment']['MESA_SHADER_CACHE_DISABLE'] == 'true'
xml = WORK / 'control/gtest.xml'
x = ET.parse(xml).getroot()
assert int(x.attrib['tests']) == 1
assert all(int(x.attrib.get(k, '0')) == 0 for k in ('failures', 'errors', 'disabled'))
assert not x.findall('.//skipped')
cases = x.findall('.//testcase')
assert len(cases) == 1
case = cases[0]
assert case.attrib['classname'] == 'VulkanBackend'
assert case.attrib['name'] == 'IndependentPairCompletesDistinctKernelsAndRejectsSharedBuffers'
assert case.attrib['result'] == 'completed'
props = {p.attrib['name']: p.attrib['value'] for p in case.findall('./properties/property')}
assert props['buffer_mutation_epochs'] == '8,12,8;8,12,8'
link = (BIN / 'link.map').read_text()
assert 'libsirius_backend.a' not in link
assert 'vulkan_smoke_test.cpp.o' in link and 'vulkan_device.cpp.o' in link
survivors = json.loads((WORK / 'survivors-before.json').read_text())
assert survivors == json.loads((WORK / 'survivors-after.json').read_text())
for name, expected in survivors.items():
    p = ROOT / name
    assert {'bytes': p.stat().st_size, 'sha256': digest(p)} == expected, name
stock = json.loads((WORK / 'protected-stock-products.json').read_text())
assert len(stock) == 33
for name, expected in stock.items():
    p = ROOT / 'bin/linux-gcc/src/sirius/backend/retained' / name
    assert {'bytes': p.stat().st_size, 'sha256': digest(p)} == expected, name
result = {
    'accepted': True, 'executed_source_revision': EXECUTED,
    'current_source_revision': subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode().strip(),
    'actions': reports, 'xml_sha256': digest(xml),
    'executable': {'bytes': (BIN / 'mutation-controls').stat().st_size, 'sha256': digest(BIN / 'mutation-controls')},
    'buffer_epochs': props['buffer_mutation_epochs'],
    'canonical_and_qualification_survivors': len(survivors), 'stock_spvs_exact': len(stock),
    'limits': 'Finite prerequisite only. No sealed owner, rendering gain, full product rebuild or scientific/release qualification. Failure/pending/overflow conclusions are source-reviewed.'
}
(WORK / 'finite-control-verification.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: result[k] for k in ('accepted', 'executed_source_revision', 'buffer_epochs', 'stock_spvs_exact', 'limits')}, indent=2))
