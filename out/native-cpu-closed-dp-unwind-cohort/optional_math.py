"""Read first shared-value PCs in the original genuine-job C8; no ABI guessing."""
from pathlib import Path
import hashlib
import json
import subprocess


def check(root,action,observation):
    def pin(path):
        with path.open('rb') as stream:
            return dict(bytes=path.stat().st_size,
                        sha256=hashlib.file_digest(stream,'sha256').hexdigest())
    path=root/'out/native-cpu-closed-dp-unwind-provider/provider/linked-native-proof.json'
    proof=json.loads(path.read_text());library=Path(action['runtime_provider'])
    assert proof['passed'] and root/proof['library']==library
    assert pin(library)==proof['library_identity']
    interface=proof['shared_value_interface']
    assert interface['passed'] and interface['core_count']==9
    assert interface['raw_core_calls']==dict(transport=2701,endpoint=3959)
    closed=proof['closed_dp_interface']
    assert closed['passed'] and closed['core_count']==9
    assert closed['raw_core_calls']==dict(dopri_phase=1576)
    assert all(item[key] is True for item in (interface,closed) for key in (
        'fixed_signatures_preserved','no_ordinary_writes','no_FP_observers','strict_scope_source_bound'))
    assert all(proof['fp_scope'][key] is True for key in (
        'saved_full_mxcsr','traps_masked','round_nearest','daz_ftz_off',
        'restored_full_mxcsr','actual_callback_code_join'))
    raw_names={'transport':'sirius_scope_transport_packet',
               'endpoint':'sirius_scope_endpoint_packet','dopri_phase':'sirius_closed_dp_packet'}
    raw_stages={name:stage for stage,name in raw_names.items()}
    selected={raw_stages[row['function']] for row in observation['simd_calls']}
    assert {'endpoint','dopri_phase'}<=selected
    intervals=observation['optional_intervals']
    assert set(intervals)==selected
    assert all(row['closed'] and row['raw_entry_sequence']<row['raw_return_sequence']
               for row in intervals.values())
    raw=dict(zip(('transport','endpoint','dopri_phase'),proof['raw_v4_objects']))
    suffixes=('dense_add8','dense_multiply8','dense_divide8','dense_sqrt8',
        'tagged_add8','tagged_subtract8','tagged_multiply8','tagged_divide8','tagged_sqrt8')
    exports={prefix+name for prefix in ('sirius_shared_value_','sirius_local_value_') for name in suffixes}
    math_objects=[proof['shared_value_math_object'],proof['closed_dp_math_object']]
    assert all(item['all_retained_function_opcodes_match'] for item in math_objects)
    math_functions={row['symbol']:row for item in math_objects for row in item['functions'] if row['symbol'] in exports}
    assert set(math_functions)==exports
    stage_prefix={'transport':'sirius_shared_value_','endpoint':'sirius_shared_value_',
                  'dopri_phase':'sirius_local_value_'}
    symbols={}
    for line in subprocess.check_output(['/usr/bin/nm','-S','-a',str(library)],text=True).splitlines():
        words=line.split()
        if len(words)==4 and (words[3] in exports or words[3]=='__cpu_features2'):
            assert words[3] not in symbols;symbols[words[3]]=int(words[0],16)
    assert set(symbols)==exports|{'__cpu_features2'}
    assert all(symbols[name]==row['linked_address'] for name,row in math_functions.items())
    rows=observation['optional_math_calls'];observed={row['stage'] for row in rows}
    assert len(rows)==len(observed) and observed==selected
    parents={row['entry_sequence']:row for row in observation['simd_parent_calls']}
    for row in rows:
        assert row['function'] in exports and row['stage']==row['expected_stage']
        assert row['function'].startswith(stage_prefix[row['stage']])
        assert row['returned'] and row['completed']
        interval=intervals[row['stage']]
        assert interval['raw_entry_sequence']<row['entry_sequence']<row['return_sequence']<interval['raw_return_sequence']
        assert row['thread']==interval['raw_thread']==row['return_thread']
        assert row['mxcsr']&0xe040==row['return_mxcsr']&0xe040==0
        assert row['mxcsr']&0x1f80==row['return_mxcsr']&0x1f80==0x1f80
        assert row['mxcsr']&~0x3f==row['return_mxcsr']&~0x3f
        assert row['library']==str(library)
        bases={int(line.split()[0].split('-')[0],16)-int(line.split()[2],16)
               for line in row['provider_maps']}
        base=row['pc']-symbols[row['function']]
        assert base==row['provider_base'] and base in bases
        assert row['cpu_features2_address']==base+symbols['__cpu_features2']
        data=bytes.fromhex(row['cpu_features2_first12']);assert len(data)==12 and data[8]&4
        assert row['raw_caller'] in row['caller_frames'] and len(row['caller_frames'])<=16
        pc=row['raw_caller']['pc']-base
        assert any(item['linked_address']<=pc<item['linked_address']+item['bytes'] or
                   item['linked_address']<=pc-1<item['linked_address']+item['bytes']
                   for item in raw[row['stage']]['functions'])
        parent=parents[row['parent_entry_sequence']]
        assert parent['entry_sequence']==interval['parent_entry_sequence']
        assert parent['entry_sequence']<row['entry_sequence']<row['return_sequence']<parent['return_sequence']
        assert parent['function']=='sirius_native_'+row['stage']+'_rows_admitted_packet'
        assert row['dispatch_job']==parent['dispatch_job']
        job=row['dispatch_job'];assert job['native_preflight_passed'] and not job['preflight']
        assert job['thread']==row['thread'] and job['stage']==('transport','endpoint','dopri_phase').index(row['stage'])+1
        assert any(frame['name']==parent['function'] for frame in row['caller_frames'])
        assert any(frame['name']==job['caller_function'] and frame['pc']==job['caller_pc'] for frame in row['caller_frames'])
        assert any(frame['name']==job['callback_function'] and frame['pc']==job['callback_pc'] for frame in row['caller_frames'])
    return dict(selected_simd_stages=sorted(selected),observed_optional_math_stages=sorted(observed),
        unobserved_optional_math_stages=sorted(selected-observed),all_selected_optional_math_observed=observed==selected,
        provider_linked_proof=pin(path),genuine_private_FP_scope_observed=True,
        shared_core_functions_observed=sorted({row['function'] for row in rows}),
        component_numerical_qualification_transferred=False,
        scope='First shared-value core entry/return in each first selected new scoped raw interval only, inside the actual same-thread positive live-job callback and statically joined full MXCSR save/mask/restore scope. Exact core/raw/provider/CPU-cache PCs and unchanged original C8/all32 assertions are separate evidence. RDI return storage and FP payloads are not read or fabricated; no inferior calls or memory/register/FP/feature writes. Intermediate sticky flags need not match outside the original quarantine. No older-hardware, benefit/all8/frame/scientific/release qualification is claimed.')
