"""Verify the first actual frame routes; never issue a performance/frame verdict."""
from pathlib import Path
import hashlib
import json
import subprocess

R=Path.cwd()
W=Path(__file__).resolve().parent
doc=lambda p:json.loads(p.read_text())


def pin(path):return {'bytes':path.stat().st_size,'sha256':hashlib.file_digest(path.open('rb'),'sha256').hexdigest()}


owner=doc(W/'observe/owner.json')
action=doc(W/'actions.json')['observe']
observation=doc(W/'observe/native-observation.json')
assert owner['status']=='terminal' and owner['passed'] and owner['child_exit']==0
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
assert all(row['absent'] for row in owner['observed_births'])
before=doc(W/'observe/inputs-before.json')
assert before==doc(W/'observe/inputs-after.json')
for name,value in before.items():assert pin(R/name)=={k:value[k] for k in ('bytes','sha256')},name
assert owner['limits']=={'seconds':45,'rss_kib':4194304,'reserve_kib':2097152}
assert observation['capture_complete'] and not observation['failures']
assert len(observation['provider_loads'])>=1 and len(observation['preflight_calls'])==3
assert len(observation['calls'])==len(observation['simd_calls'])==3
provider=Path(action['runtime_provider'])
symbols={}
for line in subprocess.check_output(['/usr/bin/nm','-S','-a',str(provider)],text=True).splitlines():
    words=line.split()
    if len(words)==4 and words[2] in ('T','t'):
        symbols.setdefault(words[3],[]).append((int(words[0],16),int(words[1],16)))
maps_by_pid={row['pid']:row['provider_maps'] for row in observation['preflight_calls']}


def entry(row):
    assert row['returned'] and Path(row['library']).resolve()==provider.resolve()
    assert row['thread']==row['return_thread'] and row['entry_sequence']<row['return_sequence']
    symbol,=symbols[row['function']]
    mappings=row.get('provider_maps',maps_by_pid[row['pid']])
    bases={int(line.split()[0].split('-')[0],16)-int(line.split()[2],16) for line in mappings}
    assert row['pc']-symbol[0] in bases
    if row['category']=='preflight_calls':assert row['return_mxcsr']==row['mxcsr']
    else:assert row['mxcsr']&0xe040==row['return_mxcsr']&0xe040==0


for row in observation['preflight_calls']+observation['calls']+observation['simd_calls']+observation['dense_calls']:
    entry(row)
for row in observation['preflight_calls']+observation['calls']:
    assert row['handled'] and row['capacity']==128
    assert 1<=row['grid_x']<=row['capacity'] and row['y']==row['z']==0
for row in observation['simd_calls']:assert row['completed'] and row['groups']==8 and row['capacity']==128
names={'sirius_native_dense_add8','sirius_native_dense_multiply8',
       'sirius_native_dense_values_add8','sirius_native_dense_values_multiply8'}
assert {row['function'] for row in observation['dense_calls']}==names
assert all(row['mask']==255 for row in observation['dense_calls'])
for row in observation['policy_calls']:
    assert row['returned'] and row['selected']==(row['workers']==0 or row['rows']>row['workers'])
result={'actual_original_frame_native_and_CPP_routes_confirmed':True,
        'owner':pin(W/'observe/owner.json'),'observation':pin(W/'observe/native-observation.json'),
        'provider':{'path':str(provider),**pin(provider)},'source':owner['source'],
        'first_preflight_arguments':[{k:row[k] for k in ('function','capacity','grid_x','groups','handled')} for row in observation['preflight_calls']],
        'actual_transport_policy':observation['policy_calls'],
        'completed_SIMD_stages':3,'active_CPP_functions':sorted(names),
        'all_first_CPP_masks':255,'qualifying_frame_completed':False,
        'runtime_speed_or_causal_share_claimed':False,'independent_math_claimed':False,
        'scope':'One deliberate first-route diagnostic of the original frame test after software warm-up. Same genuine P ELF/AD provider, original scene and cache policy. Native3/SIMD3 and all4 dense/packing functions executed and returned under IEEE32 scope. The inferior stopped after coverage; this is neither a frame retry/pass nor a timing or whole-frame causal study. It disproves an entirely unreached native/C++ path for this observed prefix only.'}
with (W/'route-result.json').open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
print(json.dumps({'native3_SIMD3_CPP4_actual_frame_routes_confirmed':True,'qualifying_frame_completed':False}))
