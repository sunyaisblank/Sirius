"""Passive first rows/pure/SIMD capture; no inferior calls or state writes.

The first actual raw SIMD entry is joined to its own concurrent parent. No
candidate return is discarded or sampled again after a refusal or exception.
"""
from pathlib import Path
import gdb
import json
import struct

work = Path.cwd() / 'out/native-cpu-upload-aware-cohort'
settings = json.loads((work / 'actions.json').read_text())['control-8']
folder = work / 'control-8'
provider = Path(settings['runtime_provider']).resolve(strict=True)
observation = {'calls': [], 'preflight_calls': [], 'simd_calls': [],
               'simd_parent_calls': [], 'capability_calls': [],
               'failures': [], 'inferior_exit': None}
armed = False
breakpoints = []
parents = {}
parent_breakpoints = {}
sequence = 0


def save():
    (folder / 'native-observation.json').write_text(json.dumps(observation, indent=2) + '\n')


def event_number():
    global sequence
    sequence += 1
    return sequence


def reg(name):
    return int(gdb.parse_and_eval('$' + name))


class Returned(gdb.FinishBreakpoint):
    def __init__(self, row):
        super().__init__(gdb.newest_frame(), internal=True)
        self.row = row

    def stop(self):
        try:
            row = self.row
            row.update(returned=True, return_sequence=event_number(),
                       return_thread=gdb.selected_thread().global_num, return_mxcsr=reg('mxcsr'))
            if row['category'] == 'simd_calls':
                row['completed'] = True  # Raw SIMD export returns void.
            elif row['category'] == 'capability_calls':
                row['supported'] = bool(reg('rax') & 0xff)
            else:
                row['handled'] = bool(reg('rax') & 0xff)
            save()
            if row['category'] not in ('simd_calls', 'capability_calls'):
                assert row['handled'], 'first observed packet refused: ' + row['function']
            assert row['thread'] == row['return_thread'], 'return thread mismatch'
            if row['category'] in ('preflight_calls', 'capability_calls'):
                assert row['return_mxcsr'] == row['mxcsr'], 'pure admission changed FP state'
            else:
                assert row['return_mxcsr'] & 0xe040 == 0, 'return IEEE32 state mismatch'
            return False
        except Exception as error:
            observation['failures'].append(self.row['function'] + ': ' + repr(error))
            save()
            return True


class CapabilityEntered(gdb.Breakpoint):
    def __init__(self):
        address = int(gdb.parse_and_eval('(void *)&sirius_native_simd_supported'))
        super().__init__('*' + hex(address), internal=True)

    def stop(self):
        self.enabled = False
        row = {'function': 'sirius_native_simd_supported',
               'category': 'capability_calls', 'entry_sequence': event_number(),
               'pc': reg('rip'), 'library': gdb.solib_name(reg('rip')),
               'pid': gdb.selected_inferior().pid,
               'thread': gdb.selected_thread().global_num,
               'mxcsr': reg('mxcsr'), 'returned': False}
        observation['capability_calls'].append(row)
        save()
        try:
            assert row['library'] and Path(row['library']).resolve(strict=True) == provider
            row['provider_maps'] = [line for line in Path(f"/proc/{row['pid']}/maps").read_text().splitlines()
                                    if str(provider) in line]
            assert row['provider_maps']
            Returned(row)
            save()
            return False
        except Exception as error:
            observation['failures'].append(row['function'] + ': ' + repr(error))
            save()
            return True


class Entered(gdb.Breakpoint):
    def __init__(self, slug, category):
        self.slug, self.category = slug, category
        self.name = ('sirius_simd_' + slug + '_packet' if category == 'simd_calls' else
                     'sirius_native_' + slug + ('_admit_packet' if category == 'preflight_calls' else '_rows_try_packet'))
        address = int(gdb.parse_and_eval('(void *)&' + self.name))
        super().__init__('*' + hex(address), internal=True)

    def stop(self):
        if self.category != 'simd_parent_calls':
            self.enabled = False
        try:
            inferior = gdb.selected_inferior()
            if self.category == 'simd_parent_calls':
                groups = struct.unpack('<I', inferior.read_memory(reg('rsp') + 40, 4).tobytes())[0]
                if groups != 8:
                    return False
            row = {'function': self.name, 'category': self.category,
                   'entry_sequence': event_number(), 'pc': reg('rip'),
                   'library': gdb.solib_name(reg('rip')), 'pid': inferior.pid,
                   'thread': gdb.selected_thread().global_num, 'input': reg('rdi'),
                   'output': reg('rdx'), 'first_x': reg('r8d') & 0xffffffff,
                   'mxcsr': reg('mxcsr'), 'returned': False}
            if self.category == 'simd_calls':
                row.update(input_words=reg('esi') & 0xffffffff,
                           output_words=reg('ecx') & 0xffffffff,
                           groups=reg('r9d') & 0xffffffff)
                caller = gdb.newest_frame().older()
                row.update(caller_function=caller.name(), caller_pc=int(caller.pc()))
            else:
                stack = struct.unpack('<5Q', inferior.read_memory(reg('rsp') + 8, 40).tobytes())
                row.update(input_bytes=reg('esi') & 0xffffffff, output_bytes=reg('ecx') & 0xffffffff,
                           y=reg('r9d') & 0xffffffff,
                           **dict(zip(('z', 'grid_x', 'grid_y', 'grid_z', 'groups'),
                                      (word & 0xffffffff for word in stack))))
            observation[self.category].append(row)
            save()  # Preserve raw arguments before all provider/range assertions.
            assert row['library'] and Path(row['library']).resolve(strict=True) == provider, 'provider mismatch'
            assert (row['input_words'] > 1 and row['output_words'] > 1 if self.category == 'simd_calls'
                    else row['input_bytes'] > 4 and row['output_bytes'] > 4), 'non-positive packet extents'
            row['capacity'] = struct.unpack('<I', inferior.read_memory(row['input'], 4).tobytes())[0]
            row['provider_maps'] = [line for line in Path(f'/proc/{inferior.pid}/maps').read_text().splitlines()
                                    if str(provider) in line]
            save()
            assert row['provider_maps'], 'actual provider mapping absent'
            if self.category != 'preflight_calls':
                assert row['mxcsr'] & 0xe040 == 0, 'entry IEEE32 state mismatch'
            if self.category == 'simd_calls':
                assert row['groups'] == 8, 'raw SIMD must execute exactly8 rows'
                parent = parents[(row['thread'], self.slug)]
                row['parent_entry_sequence'] = parent['entry_sequence']
                parent['contains_first_simd'] = True
                parent_breakpoints[self.slug].enabled = False
                for key in ('pid', 'thread', 'input', 'output', 'first_x', 'capacity', 'groups'):
                    assert row[key] == parent[key], 'raw SIMD parent mismatch: ' + key
                assert row['input_words'] * 4 == parent['input_bytes'], 'raw input word/byte mismatch'
                assert row['output_words'] * 4 == parent['output_bytes'], 'raw output word/byte mismatch'
                assert row['caller_function'] == parent['function'], 'raw caller differs from rows wrapper'
            else:
                assert 1 <= row['grid_x'] <= row['capacity'], 'active grid exceeds capacity'
                assert row['y'] == row['z'] == 0 and row['grid_y'] == row['grid_z'] == 1, 'unexpected Y/Z range'
                assert 1 <= row['groups'] <= 8, 'group count outside1..8'
                assert row['first_x'] + row['groups'] <= row['grid_x'], 'packet crosses active prefix'
                if self.category == 'simd_parent_calls':
                    parents[(row['thread'], self.slug)] = row
            Returned(row)
            save()
            return False
        except Exception as error:
            observation['failures'].append(self.name + ': ' + repr(error))
            save()
            return True


def loaded(event):
    global armed
    if armed or not event.new_objfile.filename or Path(event.new_objfile.filename).resolve() != provider:
        return
    breakpoints.append(CapabilityEntered())
    for slug in ('transport', 'endpoint', 'dopri_phase'):
        for category in ('calls', 'preflight_calls', 'simd_calls', 'simd_parent_calls'):
            breakpoint = Entered(slug, category)
            breakpoints.append(breakpoint)
            if category == 'simd_parent_calls':
                parent_breakpoints[slug] = breakpoint
    armed = True


def exited(event):
    observation['inferior_exit'] = getattr(event, 'exit_code', None)
    save()


gdb.events.new_objfile.connect(loaded)
gdb.events.exited.connect(exited)
save()
