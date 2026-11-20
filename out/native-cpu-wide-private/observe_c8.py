"""Read-only first-packet and FP-scope observations for original C8 in GDB."""

from pathlib import Path
import json
import os
import struct

import gdb


TASK = Path(os.environ['SIMD16_C8_TASK'])
OUTPUT = Path(os.environ['SIMD16_C8_OUTPUT']) / 'native-observation.json'
PROOF = json.loads((TASK / 'linked-provider.json').read_text())
PROVIDER = str(Path(PROOF['provider']['path']).resolve())
STAGES = ('transport', 'endpoint', 'dopri_phase')
WIDE_DP = 'sirius_wide_strict_dp_packet'
RECORD = {'packets': {}, 'wide_dp': {}, 'failures': [], 'exit_code': None}
BREAKPOINTS = []


def reg(name, frame=None):
    return int((frame or gdb.newest_frame()).read_register(name))


def words(address, fmt):
    return struct.unpack(fmt, gdb.selected_inferior().read_memory(address, struct.calcsize(fmt)).tobytes())


def save():
    OUTPUT.write_text(json.dumps(RECORD, indent=2) + '\n')


def ieee32():
    mxcsr = reg('mxcsr')
    assert mxcsr & 0xffc0 == 0x1f80
    return mxcsr


class Returned(gdb.FinishBreakpoint):
    def __init__(self, frame, row, entry, restores=None):
        super().__init__(frame, internal=True)
        self.thread = gdb.selected_thread().global_num
        self.row = row
        self.entry = entry
        self.restores = restores

    def stop(self):
        try:
            assert gdb.selected_thread().global_num == self.thread
            if self.restores is not None:
                assert reg('mxcsr') == self.restores
                self.row['callback_restored_mxcsr'] = self.restores
            else:
                self.row['return_mxcsr'] = ieee32()
                if self.entry:
                    assert reg('eax') & 0xff == 1
                    self.row['handled'] = True
            save()
            return False
        except Exception as error:
            RECORD['failures'].append(repr(error))
            save()
            return True


class WideDP(gdb.Breakpoint):
    def __init__(self, thread, row):
        self.row = row
        address = int(gdb.parse_and_eval('(void *)&' + WIDE_DP))
        super().__init__('*' + hex(address), internal=True)
        self.thread = thread

    def stop(self):
        self.enabled = False
        try:
            frame = gdb.newest_frame()
            caller = frame.older()
            assert caller and caller.name() == self.row['function']
            assert gdb.solib_name(frame.pc()) == PROVIDER
            assert gdb.selected_thread().global_num == self.thread
            for register, value in (
                ('rdi', self.row['input']), ('esi', self.row['input_bytes'] // 4),
                ('rdx', self.row['output']), ('ecx', self.row['output_bytes'] // 4),
                ('r8d', self.row['x']), ('r9d', 16),
            ):
                assert reg(register) == value, register
            wide = RECORD['wide_dp']
            wide.update(pc=frame.pc(), parent_pc=caller.pc(), mxcsr=ieee32(), thread=self.thread)
            BREAKPOINTS.append(Returned(frame, wide, False))
            save()
            return False
        except Exception as error:
            RECORD['failures'].append(repr(error))
            save()
            return True


class Packet(gdb.Breakpoint):
    def __init__(self, stage, wide):
        self.stage = stage
        self.name = 'sirius_native_' + stage + '_rows_' + ('wide_' if wide else '') + 'admitted_packet'
        address = int(gdb.parse_and_eval('(void *)&' + self.name))
        super().__init__('*' + hex(address), internal=True)

    def stop(self):
        try:
            frame = gdb.newest_frame()
            mode = frame.older()
            callback = mode.older()
            assert mode.name() == 'lp_native_cpu_packet_mode' and callback.name() == 'cs_exec_fn'
            assert gdb.solib_name(frame.pc()) == PROVIDER
            base = frame.pc() - PROOF['packet_entries'][self.name]['entry']
            assert mode.pc() == base + PROOF['packet_entries'][self.name]['return']
            assert callback.pc() == base + PROOF['callback']['mode_return']
            job = reg('rbx', callback)
            grid, = words(job, '<I')
            admitted, width, stage = words(job + 63, '<BBB')
            x = reg('r8d')
            z, grid_x, gy, gz, groups = words(reg('rsp') + 8, '<I4xI4xI4xI4xI')
            assert admitted == 1 and stage == STAGES.index(self.stage) + 1
            assert grid == grid_x and x + groups <= grid and (reg('r9d'), z, gy, gz) == (0, 0, 1, 1)
            current, = words(job + 72, '<Q')
            table, elements = words(current + 24, '<QI')
            assert elements == 32
            input_pointer, input_bytes = words(table, '<QI')
            output_pointer, output_bytes = words(table + 64, '<QI')
            assert (input_pointer, input_bytes, output_pointer, output_bytes) == (
                reg('rdi'), reg('esi'), reg('rdx'), reg('ecx'))
            key = f'{self.stage}:{grid}:{groups}'
            if key in RECORD['packets']:
                return False
            thread = gdb.selected_thread().global_num
            row = dict(function=self.name, pc=frame.pc(), mode_pc=mode.pc(), callback_pc=callback.pc(),
                       job=job, admitted=admitted, width=width, stage=stage, grid=grid, x=x, groups=groups,
                       input=input_pointer, input_bytes=input_bytes, output=output_pointer,
                       output_bytes=output_bytes, capacity=words(input_pointer, '<I')[0],
                       thread=thread, mxcsr=ieee32())
            RECORD['packets'][key] = row
            BREAKPOINTS.append(Returned(frame, row, True))
            if groups == 16 and self.stage == 'dopri_phase':
                saved_mxcsr, = words(reg('rsp', callback) + 4, '<I')
                row['callback_saved_mxcsr'] = saved_mxcsr
                BREAKPOINTS.extend((Returned(callback, row, False, saved_mxcsr), WideDP(thread, row)))
            save()
            return False
        except Exception as error:
            RECORD['failures'].append(repr(error))
            save()
            return True


def loaded(event):
    if event.new_objfile.filename != PROVIDER:
        return
    for stage in STAGES:
        for wide in (False, True):
            BREAKPOINTS.append(Packet(stage, wide))


def exited(event):
    RECORD['exit_code'] = getattr(event, 'exit_code', None)
    save()


gdb.events.new_objfile.connect(loaded)
gdb.events.exited.connect(exited)
