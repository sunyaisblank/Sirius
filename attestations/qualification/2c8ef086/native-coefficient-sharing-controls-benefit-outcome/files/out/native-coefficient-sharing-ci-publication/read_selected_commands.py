import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
STAGE = ROOT / 'bin/windows-msvc/native-coefficient-sharing-matched-current'
audit = json.loads((ROOT / 'out/native-coefficient-sharing-pe-audit/actual-e1c618f8-audit.json').read_text())


def text(path):
    data = path.read_bytes()
    return data.decode('utf-16' if data[:2] in (b'\xff\xfe', b'\xfe\xff') else 'utf-8-sig')


def selected(paths, suffix):
    records = []
    for path in paths:
        keep = False
        for line in text(path).splitlines():
            if line.startswith('^'):
                keep = line.replace('\\', '/').lower().endswith(suffix.lower())
            elif keep and line.strip():
                records.append(line)
    return records


def normalize(command):
    tokens = re.findall(r'(?:[^\s"]|"[^"]*")+', command)
    return [word.replace('\\', '/').lower() for word in tokens
            if not word.lower().startswith(('/fo', '/fd'))
            and 'out\\native-coefficient-sharing-producer-inputs\\' not in word.lower()]


backend = STAGE / 'metadata/canonical-backend/src/sirius/backend/sirius_backend.dir/Release'
canonical = selected(list(backend.rglob('CL.command.*.tlog')), 'src/sirius/backend/retained_compute.cpp')
assert len(canonical) == 1
result = {'canonical_compute_command': canonical[0], 'private': {}, 'canonical_Vulkan_MAP_functions': {}}
for target, row in audit['private_compiles'].items():
    lines = [line for name in row['command_tlogs'] for line in text(STAGE / name).splitlines()
             if line.strip() and not line.startswith('^')]
    assert len(lines) == 1
    command = lines[0]
    for flag in ('/W4', '/WX', '/O2', '/Ob2', '/MD', '/fp:precise', '/std:c++latest',
                 '/permissive-', 'SIRIUS_CONTRACT_MODE=2'):
        assert flag.lower() in command.lower()
    if target.endswith('_compute'):
        assert normalize(canonical[0]) == normalize(command)
    result['private'][target] = {'command': command,
        'canonical_compute_usage_exact_after_only_private_I_Fo_Fd_removal': target.endswith('_compute')}
links = {}
for arm, row in audit['products'].items():
    assert len(row['link_command_tlogs']) == 1
    # A LINK command has continuation object-path lines in its tlog record.
    command = '\n'.join(line for line in text(STAGE / row['link_command_tlogs'][0]).splitlines()
                        if line.strip() and not line.startswith('^'))
    assert '/verbose:lib' in command.lower() and '/wholearchive' not in command.lower()
    links[arm] = command
    mapping = text(STAGE / ('payload/' + arm + '/original_' + arm + '_controls.map'))
    functions = {}
    for name in ('LoadKernel', 'GetOrCreatePipeline', 'CreateBuffer', 'DispatchCommands'):
        lines = [line for line in mapping.splitlines()
                 if re.search(r'\s\?' + name + r'@VulkanDevice@backend@sirius@@', line)]
        assert len(lines) == 1 and 'sirius_backend:vulkan_device.obj' in lines[0]
        functions[name] = lines[0]
    result['canonical_Vulkan_MAP_functions'][arm] = functions
assert links['baseline'].replace('baseline', 'candidate').replace('BASELINE', 'CANDIDATE') == links['candidate']
result.update(passed=True, link_commands=links,
    scope='Actual5CL+2LINK and exact canonical compute usage/4Vulkan MAP functions per arm. Mature unique archive object and LINK.read checks separately required; no individual Loaded extraction line asserted.',
    passive_prewrite_correction='Initial inline reader expected one physical LINK line; actual single LINK tlog record has two continuation object lines. Join complete record. No build or workload repeated.')
output = WORK / 'root-selected-command-details.json'
assert not output.exists()
output.write_text(json.dumps(result, indent=2) + '\n')
print('Actual5CL/2LINK strict; both compute usage lists exact canonical; four Vulkan functions per arm select canonical archive.')
