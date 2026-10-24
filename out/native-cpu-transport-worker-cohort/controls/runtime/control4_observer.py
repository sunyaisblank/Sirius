"""Observe the next original C4 once, including one fixed-time stack sample.

No inferior calls, register/data writes, altered test inputs or retry. The
debugger remains inside the original 300s owner. Interrupt/resume overhead is
included; this numerical control supplies no performance-benefit evidence.
The owning Python process sends one SIGINT to its own debugger child. A small
independent sleep-process check rejected an internal GDB Python timer because
its interrupt was deferred, and verified external interrupt/capture/resume.
No extra Sirius test is used to check the observation harness.
"""
from pathlib import Path
import gdb
import json
import struct
import time

work = Path.cwd() / 'out/native-cpu-transport-worker-cohort/controls/runtime'
folder = work / 'control-4'
action = json.loads((work / 'actions.json').read_text())['control-4']
provider = Path(action['runtime_provider']).resolve(strict=True)
started = time.monotonic()
record = {'calls': [], 'failures': [], 'sample_after_seconds': 45,
          'sample': None, 'inferior_exit': None, 'last_stop_signal': None,
          'scope': 'One original C4 under bounded observation; all overhead included. No performance evidence.'}
armed = False
breakpoints = []


def save():
    (folder / 'route-observation.json').write_text(json.dumps(record, indent=2) + '\n')


def reg(name):
    return int(gdb.parse_and_eval('$' + name))


def frames(limit=8):
    result, frame = [], gdb.newest_frame()
    while frame is not None and len(result) < limit:
        pc = int(frame.pc())
        result.append({'name': frame.name(), 'pc': pc, 'library': gdb.solib_name(pc)})
        frame = frame.older()
    return result


class Returned(gdb.FinishBreakpoint):
    def __init__(self, row):
        super().__init__(gdb.newest_frame(), internal=True)
        self.row = row

    def stop(self):
        self.row.update(returned=True, return_seconds=time.monotonic() - started,
                        return_thread=gdb.selected_thread().global_num)
        if self.row['category'] == 'native_bool':
            self.row['handled'] = bool(reg('rax') & 255)
        save()
        return False

    def out_of_scope(self):
        self.row['finish_out_of_scope'] = True
        save()


class Entered(gdb.Breakpoint):
    def __init__(self, name, category, maximum=1):
        self.name, self.category, self.maximum, self.count = name, category, maximum, 0
        address = int(gdb.parse_and_eval('(void *)&' + name))
        super().__init__('*' + hex(address), internal=True)

    def stop(self):
        try:
            self.count += 1
            if self.count == self.maximum:
                self.enabled = False
            row = {'function': self.name, 'category': self.category,
                   'entry_seconds': time.monotonic() - started, 'returned': False,
                   'thread': gdb.selected_thread().global_num, 'pc': reg('rip'),
                   'library': gdb.solib_name(reg('rip')), 'frames': frames()}
            record['calls'].append(row)
            assert row['library'] and Path(row['library']).resolve(strict=True) == provider
            inferior = gdb.selected_inferior()
            if self.category == 'native_bool':
                stack = struct.unpack('<5Q', inferior.read_memory(reg('rsp') + 8, 40).tobytes())
                row.update(first_x=reg('r8d') & 0xffffffff,
                           **dict(zip(('z', 'grid_x', 'grid_y', 'grid_z', 'groups'),
                                      (value & 0xffffffff for value in stack))))
                row['capacity'] = struct.unpack('<I', inferior.read_memory(reg('rdi'), 4).tobytes())[0]
            elif self.category == 'raw_SIMD':
                row['groups'] = reg('r9d') & 0xffffffff
                row['unexpected_for_original_single_row_C4'] = True
            save()
            Returned(row)
            return False
        except Exception as error:
            record['failures'].append(self.name + ': ' + repr(error))
            save()
            return True


def loaded(event):
    global armed
    if armed or not event.new_objfile.filename:
        return
    if Path(event.new_objfile.filename).resolve() != provider:
        return
    armed = True
    for slug in ('transport', 'endpoint', 'dopri_phase'):
        for suffix in ('_admit_packet', '_rows_try_packet', '_try_packet'):
            breakpoints.append(Entered('sirius_native_' + slug + suffix, 'native_bool'))
        breakpoints.append(Entered('sirius_simd_' + slug + '_packet', 'raw_SIMD'))
    for name in ('lvp_shader_compile', 'llvmpipe_update_cs_variant', 'gallivm_compile_module'):
        breakpoints.append(Entered(name, 'ordinary_compilation', maximum=6))
    save()


def exited(event):
    record['inferior_exit'] = getattr(event, 'exit_code', None)
    record['exit_seconds'] = time.monotonic() - started
    save()


def resume_after_sample():
    if record['inferior_exit'] is not None:
        return
    if not (folder / 'sample-request.json').is_file() or record['last_stop_signal'] != 'SIGINT':
        record['failures'].append('Unexpected inferior stop before fixed-time sample')
        save()
        return
    record['sample'] = {'elapsed_seconds': time.monotonic() - started,
                        'all_threads_bt8': gdb.execute('thread apply all bt 8', to_string=True)}
    save()
    gdb.execute('continue')


def stopped(event):
    if isinstance(event, gdb.SignalEvent):
        record['last_stop_signal'] = event.stop_signal
        save()


gdb.events.new_objfile.connect(loaded)
gdb.events.exited.connect(exited)
gdb.events.stop.connect(stopped)
save()
