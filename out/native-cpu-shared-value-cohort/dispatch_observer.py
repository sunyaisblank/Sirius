# Runs in the original passive observer namespace; no inferior calls/writes.
_dispatch_contract = json.loads((work.parent / 'native-cpu-shared-value-provider/provider/linked-native-proof.json').read_text())
_dispatch_layout = _dispatch_contract['job_observation_layout']
_dispatch_register = _dispatch_contract['callback_live_job_register']
assert _dispatch_contract['passed'] and _dispatch_register['no_intervening_job_register_writes']
assert _dispatch_register['frame'] == 'cs_exec_fn' and _dispatch_register['register'] == 'rbx'


def capture_dispatch_job(row, slug):
    mode = gdb.newest_frame().older()
    caller = mode.older() if mode else None
    assert mode and mode.name() == 'lp_native_cpu_packet_mode', 'actual Mesa mode caller unavailable'
    assert caller and caller.name() == _dispatch_register['frame'], 'actual Mesa callback unavailable'
    job_pointer = int(caller.read_register(_dispatch_register['register']))
    inferior = gdb.selected_inferior()
    raw = inferior.read_memory(job_pointer, _dispatch_layout['job_bytes']).tobytes()
    current = struct.unpack_from('<Q', raw, _dispatch_layout['current_offset'])[0]
    table = inferior.read_memory(current + _dispatch_layout['table_view_offset'], _dispatch_layout['buffer_bytes']).tobytes()
    descriptor_pointer = struct.unpack_from('<Q', table, _dispatch_layout['buffer_pointer_offset'])[0]
    table_elements = struct.unpack_from('<I', table, _dispatch_layout['buffer_extent_offset'])[0]
    assert table_elements * 4 == 2 * _dispatch_layout['descriptor_bytes']
    descriptors = inferior.read_memory(descriptor_pointer, 2 * _dispatch_layout['descriptor_bytes']).tobytes()
    def view(offset):
        return (struct.unpack_from('<Q', descriptors, offset + _dispatch_layout['buffer_pointer_offset'])[0],
                struct.unpack_from('<I', descriptors, offset + _dispatch_layout['buffer_extent_offset'])[0])
    input_pointer,input_bytes=view(0)
    output_pointer,output_bytes=view(_dispatch_layout['descriptor_bytes'])
    value = dict(job=job_pointer, native_preflight_passed=bool(raw[_dispatch_layout['native_preflight_passed_offset']]),
        preflight=False, preflight_evidence='Exact final callback zeroes fifth mode argument before the observed call',
        stage=raw[_dispatch_layout['stage_offset']],
        **dict(zip(('grid_x','grid_y','grid_z'),struct.unpack_from('<3I',raw,_dispatch_layout['grid_size_offset']))),
        input=input_pointer,input_bytes=input_bytes,output=output_pointer,output_bytes=output_bytes,
        caller_function=mode.name(), caller_pc=int(mode.pc()), callback_function=caller.name(),callback_pc=int(caller.pc()),
        thread=gdb.selected_thread().global_num, job_type_bytes=_dispatch_layout['job_bytes'],
        descriptor_type_bytes=_dispatch_layout['descriptor_bytes'],layout=_dispatch_layout,
        register_evidence=_dispatch_register,
        source='Genuine unwound callee-saved callback job at existing first private rows entry; every layout read is compile-time asserted')
    actual_callback=int(gdb.parse_and_eval('(void *)&cs_exec_fn'))
    assert Path(gdb.solib_name(actual_callback)).resolve(strict=True)==provider
    base=actual_callback-_dispatch_register['function_address']
    expected_call,=[entry for entry in _dispatch_contract['actual_hidden_admitted_worker_calls'] if entry['target']==row['function']]
    assert value['callback_pc']==base+_dispatch_register['return_address']
    assert value['caller_pc']==base+expected_call['return_address']
    assert value['native_preflight_passed'] and not value['preflight']
    assert value['stage'] == ('transport', 'endpoint', 'dopri_phase').index(slug) + 1
    for key in ('thread','grid_x','grid_y','grid_z','input','input_bytes','output','output_bytes'):
        assert value[key] == row[key], 'actual dispatch job/worker mismatch: ' + key
    row['dispatch_job'] = value
