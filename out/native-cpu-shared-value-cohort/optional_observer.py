# Inserted into the original passive observer namespace before its loaded hook.
# The value ABI's RDI is uninitialized compiler return storage, not an input bank.
import subprocess

optional_breakpoints=[]
optional_armed_stage=None
observation['optional_math_calls']=[]
observation['optional_intervals']={}
proof=json.loads((work.parent/'native-cpu-shared-value-provider/provider/linked-native-proof.json').read_text())
assert proof['passed'] and proof['shared_value_interface']['passed']
assert all(proof['fp_scope'][key] is True for key in (
    'saved_full_mxcsr','traps_masked','round_nearest','daz_ftz_off',
    'restored_full_mxcsr','actual_callback_code_join'))
scope_raw_names={'transport':'sirius_scope_transport_packet',
                 'endpoint':'sirius_scope_endpoint_packet',
                 'dopri_phase':'sirius_scope_dp_packet'}
scope_raw_stages={name:stage for stage,name in scope_raw_names.items()}
raw_ranges={}
for slug,item in zip(('transport','endpoint','dopri_phase'),proof['raw_v4_objects']):
    assert item['all_retained_function_opcodes_match']
    assert scope_raw_names[slug] in {row['symbol'] for row in item['functions']}
    raw_ranges[slug]=[(row['linked_address'],row['linked_address']+row['bytes'])
                     for row in item['functions']]
math=proof['shared_value_math_object']
assert math['all_retained_function_opcodes_match']
shared_core_names=tuple('sirius_shared_value_'+name for name in (
    'dense_add8','dense_multiply8','dense_divide8','dense_sqrt8',
    'tagged_add8','tagged_subtract8','tagged_multiply8','tagged_divide8','tagged_sqrt8'))
math_functions={row['symbol']:row for row in math['functions']
                if row['symbol'].startswith('sirius_shared_value_')}
assert set(math_functions)==set(shared_core_names)
optional_symbols={}
for line in subprocess.check_output(['/usr/bin/nm','-S','-a',str(provider)],text=True).splitlines():
    words=line.split()
    if len(words)==4 and (words[3] in shared_core_names or words[3]=='__cpu_features2'):
        assert words[3] not in optional_symbols
        optional_symbols[words[3]]=int(words[0],16)
assert set(optional_symbols)==set(shared_core_names)|{'__cpu_features2'}
assert all(optional_symbols[name]==row['linked_address'] for name,row in math_functions.items())


def arm_optional_math(row):
    global optional_armed_stage
    slug=scope_raw_stages[row['function']]
    assert slug not in observation['optional_intervals']
    parent,=[item for item in observation['simd_parent_calls']
             if item['entry_sequence']==row['parent_entry_sequence']]
    assert parent['dispatch_job']['native_preflight_passed'] is True
    observation['optional_intervals'][slug]=dict(raw_entry_sequence=row['entry_sequence'],
        raw_thread=row['thread'],parent_entry_sequence=parent['entry_sequence'],closed=False)
    optional_armed_stage=slug
    for breakpoint in optional_breakpoints:
        breakpoint.thread=row['thread'];breakpoint.enabled=True


def close_optional_math(row):
    global optional_armed_stage
    slug=scope_raw_stages[row['function']]
    interval=observation['optional_intervals'][slug]
    interval.update(closed=True,raw_return_sequence=row['return_sequence'])
    if optional_armed_stage==slug:
        for breakpoint in optional_breakpoints:breakpoint.enabled=False
        optional_armed_stage=None


class OptionalReturned(gdb.FinishBreakpoint):
    def __init__(self,row):
        super().__init__(gdb.newest_frame(),internal=True);self.row=row

    def stop(self):
        try:
            row=self.row
            row.update(returned=True,completed=True,return_sequence=event_number(),
                return_thread=gdb.selected_thread().global_num,return_mxcsr=reg('mxcsr'))
            save()
            assert row['return_thread']==row['thread']
            assert row['return_mxcsr']&0xe040==0
            assert row['return_mxcsr']&0x1f80==0x1f80
            assert (row['return_mxcsr']&~0x3f)==(row['mxcsr']&~0x3f)
            return False
        except Exception as error:
            observation['failures'].append('shared-value return: '+repr(error));save();return True


class OptionalEntered(gdb.Breakpoint):
    def __init__(self,name):
        self.name=name
        address=int(gdb.parse_and_eval('(void *)&'+name))
        super().__init__('*'+hex(address),internal=True);self.enabled=False

    def stop(self):
        global optional_armed_stage
        for breakpoint in optional_breakpoints:breakpoint.enabled=False
        expected_stage=optional_armed_stage;optional_armed_stage=None
        row=dict(function=self.name,entry_sequence=event_number(),pc=reg('rip'),
            library=gdb.solib_name(reg('rip')),pid=gdb.selected_inferior().pid,
            thread=gdb.selected_thread().global_num,mxcsr=reg('mxcsr'),
            returned=False,expected_stage=expected_stage)
        observation['optional_math_calls'].append(row);save()
        try:
            assert row['library'] and Path(row['library']).resolve(strict=True)==provider
            row['provider_maps']=[line for line in Path(f"/proc/{row['pid']}/maps").read_text().splitlines()
                                  if str(provider) in line]
            bases={int(line.split()[0].split('-')[0],16)-int(line.split()[2],16)
                   for line in row['provider_maps']}
            base=row['pc']-optional_symbols[self.name]
            assert base in bases
            row['provider_base']=base
            cache=base+optional_symbols['__cpu_features2']
            row['cpu_features2_address']=cache
            row['cpu_features2_first12']=bytes(gdb.selected_inferior().read_memory(cache,12)).hex()
            frames=[];caller=gdb.newest_frame().older()
            for _ in range(16):
                if caller is None:break
                frames.append(dict(name=caller.name(),pc=int(caller.pc())));caller=caller.older()
            row['caller_frames']=frames
            matches=[(slug,frame) for frame in frames for slug,ranges in raw_ranges.items()
                     if any(lo<=frame['pc']-base<hi or lo<=frame['pc']-base-1<hi for lo,hi in ranges)]
            assert {slug for slug,_ in matches}=={expected_stage}
            row['stage']=expected_stage;row['raw_caller']=matches[0][1]
            interval=observation['optional_intervals'][expected_stage]
            assert row['thread']==interval['raw_thread']
            parent,=[item for item in observation['simd_parent_calls']
                     if item['entry_sequence']==interval['parent_entry_sequence']]
            job=parent['dispatch_job']
            assert any(frame['name']==parent['function'] for frame in frames)
            assert any(frame['name']==job['caller_function'] and frame['pc']==job['caller_pc'] for frame in frames)
            assert any(frame['name']==job['callback_function'] and frame['pc']==job['callback_pc'] for frame in frames)
            row['parent_entry_sequence']=parent['entry_sequence'];row['dispatch_job']=dict(job)
            assert job['native_preflight_passed'] and job['thread']==row['thread']
            save()
            assert row['mxcsr']&0xe040==0 and row['mxcsr']&0x1f80==0x1f80
            assert bytes.fromhex(row['cpu_features2_first12'])[8]&4
            OptionalReturned(row)
            return False
        except Exception as error:
            observation['failures'].append('shared-value entry: '+repr(error));save();return True
