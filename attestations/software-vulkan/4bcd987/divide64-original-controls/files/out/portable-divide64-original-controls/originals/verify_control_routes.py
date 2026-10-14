"""Join Factored Endpoint route loads, successful calls and complete readbacks."""
from pathlib import Path
import hashlib, json, re, shutil, xml.etree.ElementTree as ET

R = Path.cwd().resolve()
W = Path(__file__).resolve().parent
folder = W / 'control-1'
owner = json.loads((folder / 'owner.json').read_text())
assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
props = {e.get('name'): e.get('value') for e in ET.parse(folder / 'gtest.xml').getroot().findall('.//property')}
arrays = json.loads((W / 'producer-acceptance.json').read_text())['actual_ELF_arrays']['candidate']
lines = (folder / 'stderr.log').read_text().splitlines()
blocks = []
for line in lines:
    match = re.fullmatch(r'\[FactoredEndpoint\] create: portable=(\d);wide=(\d);decline_rte32=(\d);fma_enabled=(\d)', line)
    if match:
        blocks.append({'portable': int(match[1]), 'wide': int(match[2]), 'decline': int(match[3]),
                       'fma': int(match[4]), 'lines': []})
    elif blocks:
        blocks[-1]['lines'].append(line)
assert len(blocks) == 4
assert {(b['decline'], b['wide']) for b in blocks} == {(d,w) for d in [0,1] for w in [0,1]}
expected_rows = [29,29,29,4,8,29,27,27,23,23,23,4,1,29]
records = []
for block in blocks:
    assert block['portable'] == 1 and block['fma'] == 0
    assert '[FactoredEndpoint] create returned: success=1' in block['lines']
    wide, decline = block['wide'], block['decline']
    stages = {}
    for line in block['lines']:
        match = re.fullmatch(r'\[FactoredEndpoint\] selected_stage=(\w+);sha256=([0-9a-f]{64})', line)
        if match:
            assert match[1] not in stages
            stages[match[1]] = match[2]
    assert set(stages) == {'camera','transport','endpoint','dense','initialize','ray_camera','dopri_phase'}
    suffix = 'PortableFp64Shader' if wide else 'PortableShader'
    names = {'camera': 'kCamera'+suffix, 'dense': 'kDense'+suffix, 'initialize': 'kInitialize'+suffix,
             'ray_camera': 'kRayCamera'+suffix,
             'transport': 'kTransport'+suffix if decline else 'kTransportPortableNormalSumShader',
             'endpoint': 'kEndpoint'+suffix if decline else 'kEndpointPortableNormalSumShader',
             'dopri_phase': 'kDopriPhase'+suffix if decline else 'kDopriPhasePortableNormalSumShader'}
    for stage, name in names.items():
        assert stages[stage] == arrays[name]['sha256'], (decline,wide,stage)
    selected = [line.split('=')[1] for line in block['lines'] if line.startswith('[FactoredEndpoint] selected_endpoint_sha256=')]
    assert selected == [stages['endpoint']]
    begins = []
    returns = []
    for line in block['lines']:
        m = re.fullmatch(r'\[FactoredEndpoint\] begin: call=(\d+);rows=(\d+);wide=(\d);decline_rte32=(\d)', line)
        if m:
            assert int(m[3]) == wide and int(m[4]) == decline
            begins.append((int(m[1]),int(m[2])))
        m = re.match(r'\[FactoredEndpoint\] returned: call=(\d+);success=(\d);', line)
        if m:
            assert m[2] == '1'
            returns.append(int(m[1]))
    assert begins == list(enumerate(expected_rows)) and returns == list(range(14))
    readbacks = {}
    for call, rows in begins:
        key = f'factored_endpoint_readback_{decline}_{wide}_{call}'
        value = dict(part.split('=',1) for part in props[key].split(';'))
        assert set(value) == {'rows','bytes','sha256'}
        assert int(value['rows']) == rows and int(value['bytes']) == rows*3829*4
        assert re.fullmatch('[0-9a-f]{64}', value['sha256'])
        readbacks[key] = value
    records.append({'declined_RTE32': bool(decline), 'wide': bool(wide), 'selected_endpoint_array': names['endpoint'],
                    'actual_loaded_endpoint': arrays[names['endpoint']], 'successful_calls': begins,
                    'complete_readbacks': readbacks,
                    'binding_scope': 'Actual LoadKernel word hashing and completed original Endpoint calls on source-bound same instance; no separately logged per-command handle.'})
original = R / 'attestations/software-vulkan/8a27321/vulkan-descriptor-reuse/files/out/vulkan-descriptor-reuse/control-4/gtest.xml'
oldprops = {e.get('name'):e.get('value') for e in ET.parse(original).getroot().findall('.//property')}
keys = {key for record in records for key in record['complete_readbacks']}
assert len(keys) == 56
assert {k:v for k,v in props.items() if k.startswith('factored_endpoint_readback_')} == {k:oldprops[k] for k in keys}
shutil.copyfile(original, W / 'original-factored-readbacks.xml')
result = {'passed': True, 'routes': records, 'complete_readbacks': 56,
          'complete_readbacks_exact_to_prior_original_control': True,
          'scope': 'Original finite Factored Endpoint control exercises candidate normal and default core routes, plus candidate masked-wide core route. Original reference/admission/refusal/recovery assertions passed. This does not alter actual RTE32 provider behavior or replace optional owner/native execution qualification.'}
(W / 'endpoint-route-acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'passed':True,'routes':4,'complete_readbacks':56,'candidate_core_source_bound_route':True}))
