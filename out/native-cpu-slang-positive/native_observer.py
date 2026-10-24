"""GDB passive first-success observation; sourced by the single original C8.

No inferior calls, register writes, data writes, extra dispatches, or retries.
Only the first call per stage is observed, so a refused first packet remains
explicit evidence and fails positive admission rather than sampling for success.
"""
from pathlib import Path
import gdb
import json
import struct

work = Path.cwd() / 'out/native-cpu-slang-positive'
settings = json.loads((work / 'actions.json').read_text())['control-8']
folder = work / 'control-8'
provider = Path(settings['runtime_provider']).resolve(strict=True)
observation = {'calls': [], 'failures': [], 'inferior_exit': None}
armed = False
breakpoints = []


def save():
    (folder / 'native-observation.json').write_text(json.dumps(observation, indent=2) + '\n')


def reg(name):
    return int(gdb.parse_and_eval('$' + name))


class Returned(gdb.FinishBreakpoint):
    def __init__(self, row):
        super().__init__(gdb.newest_frame(), internal=True)
        self.row = row

    def stop(self):
        try:
            self.row.update(returned=True, handled=bool(reg('rax') & 0xff),
                            return_thread=gdb.selected_thread().global_num,
                            return_mxcsr=reg('mxcsr'))
            assert self.row['handled'], 'First original packet was refused: ' + self.row['function']
            assert self.row['thread'] == self.row['return_thread']
            assert self.row['return_mxcsr'] & 0xe040 == 0
            save()
            return False
        except Exception as error:
            observation['failures'].append(repr(error))
            save()
            return True


class Entered(gdb.Breakpoint):
    def __init__(self, name):
        self.name = name
        address = int(gdb.parse_and_eval('(void *)&' + name))
        super().__init__('*' + hex(address), internal=True)

    def stop(self):
        self.enabled = False
        try:
            library = gdb.solib_name(reg('rip'))
            assert library and Path(library).resolve(strict=True) == provider
            inferior = gdb.selected_inferior()
            stack = struct.unpack('<5Q', inferior.read_memory(reg('rsp') + 8, 40).tobytes())
            row = {'function': self.name, 'pc': reg('rip'), 'library': library,
                   'pid': inferior.pid, 'thread': gdb.selected_thread().global_num,
                   'input': reg('rdi'), 'input_bytes': reg('esi') & 0xffffffff,
                   'output': reg('rdx'), 'output_bytes': reg('ecx') & 0xffffffff,
                   'first_x': reg('r8d') & 0xffffffff, 'y': reg('r9d') & 0xffffffff,
                   **dict(zip(('z', 'grid_x', 'grid_y', 'grid_z', 'groups'),
                              (word & 0xffffffff for word in stack))),
                   'mxcsr': reg('mxcsr'), 'returned': False}
            assert row['input_bytes'] > 4 and row['output_bytes'] > 4
            assert row['mxcsr'] & 0xe040 == 0
            row['rows'] = struct.unpack('<I', inferior.read_memory(row['input'], 4).tobytes())[0]
            assert row['rows'] == row['grid_x'] and row['y'] == row['z'] == 0
            assert row['grid_y'] == row['grid_z'] == 1
            assert 1 <= row['groups'] <= 8 and row['first_x'] + row['groups'] <= row['rows']
            row['provider_maps'] = [line for line in Path(f'/proc/{inferior.pid}/maps').read_text().splitlines()
                                    if str(provider) in line]
            assert row['provider_maps']
            observation['calls'].append(row)
            Returned(row)
            save()
            return False
        except Exception as error:
            observation['failures'].append(repr(error))
            save()
            return True


def loaded(event):
    global armed
    if armed or not event.new_objfile.filename:
        return
    if Path(event.new_objfile.filename).resolve() != provider:
        return
    for slug in ('transport', 'endpoint', 'dopri_phase'):
        breakpoints.append(Entered('sirius_native_' + slug + '_try_packet'))
    armed = True


def exited(event):
    observation['inferior_exit'] = getattr(event, 'exit_code', None)
    save()


gdb.events.new_objfile.connect(loaded)
gdb.events.exited.connect(exited)
save()
