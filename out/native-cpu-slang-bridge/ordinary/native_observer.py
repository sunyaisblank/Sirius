"""GDB passive native-entry/return observer for the original four-byte guard."""
# Sourced after original capture_jit/features/properties; never calls inferior code.
native_rows = []
native_armed = False


def native_write():
    (folder / 'native-calls.json').write_text(json.dumps(native_rows, indent=2) + '\n')


class NativeReturned(gdb.FinishBreakpoint):
    def __init__(self, row):
        super().__init__(gdb.newest_frame(), internal=True)
        self.row = row

    def stop(self):
        try:
            self.row['returned'] = True
            self.row['handled'] = bool(reg('rax') & 0xff)
            self.row['output_after'] = memory(self.row['output'], self.row['output_bytes']).hex()
            assert not self.row['handled']
            assert self.row['output_before'] == self.row['output_after']
            native_write()
            return False
        except Exception as error:
            return failure(error)


class NativeEntry(gdb.Breakpoint):
    def __init__(self, name):
        self.native_name = name
        address = int(gdb.parse_and_eval('(void *)&' + name))
        super().__init__('*' + hex(address), internal=True)

    def stop(self):
        try:
            actual_library = gdb.solib_name(reg('rip'))
            assert actual_library and pathlib.Path(actual_library).resolve() == pathlib.Path(os.environ['SIRIUS_ACTUAL_PROVIDER']).resolve()
            stack_args = struct.unpack('<5Q', memory(reg('rsp') + 8, 40))
            row = {'function': self.native_name, 'pc': reg('rip'), 'library': actual_library,
                   'thread': gdb.selected_thread().global_num,
                   'input': reg('rdi'), 'input_bytes': reg('esi') & 0xffffffff,
                   'output': reg('rdx'), 'output_bytes': reg('ecx') & 0xffffffff,
                   'first_x': reg('r8d') & 0xffffffff, 'y': reg('r9d') & 0xffffffff,
                   **dict(zip(('z', 'grid_x', 'grid_y', 'grid_z', 'groups'),
                              (word & 0xffffffff for word in stack_args))),
                   'mxcsr': reg('mxcsr'), 'returned': False}
            assert row['input_bytes'] == 4 and 0 < row['output_bytes'] <= 4096
            assert row['mxcsr'] & 0xe040 == 0
            row['output_before'] = memory(row['output'], row['output_bytes']).hex()
            native_rows.append(row)
            NativeReturned(row)
            native_write()
            return False
        except Exception as error:
            return failure(error)


def native_arm(event):
    global native_armed
    if native_armed or not event.new_objfile.filename:
        return
    if pathlib.Path(event.new_objfile.filename).resolve() != pathlib.Path(os.environ['SIRIUS_ACTUAL_PROVIDER']).resolve():
        return
    for slug in ('transport', 'endpoint', 'dopri_phase'):
        NativeEntry('sirius_native_' + slug + '_try_packet')
    native_armed = True


gdb.events.new_objfile.connect(native_arm)
