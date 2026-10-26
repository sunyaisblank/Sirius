# Inserted into the unchanged AN observer namespace before its loaded hook.
import subprocess

optional_breakpoints = []
optional_armed_stage = None
observation['optional_math_calls'] = []
observation['optional_intervals'] = {}
proof = json.loads((Path.cwd() / 'out/native-cpu-v4-provider/provider/linked-native-proof.json').read_text())
assert proof['passed']
raw_ranges = {}
for index, slug in enumerate(('transport', 'endpoint', 'dopri_phase'), 4):
    item = proof['objects'][index]
    assert item['object']['path'].endswith('/' + slug + '.o')
    raw_ranges[slug] = [(row['linked_address'], row['linked_address'] + row['bytes'])
                        for row in item['functions']]
optional_symbols = {}
for line in subprocess.check_output(['/usr/bin/nm', '-S', '-a', str(provider)], text=True).splitlines():
    words = line.split()
    if len(words) == 4 and (words[3].startswith('sirius_native_v4_') or words[3] == '__cpu_features2'):
        assert words[3] not in optional_symbols
        optional_symbols[words[3]] = int(words[0], 16)
assert '__cpu_features2' in optional_symbols


def arm_optional_math(row):
    global optional_armed_stage
    slug = row['function'][len('sirius_simd_'):-len('_packet')]
    assert slug not in observation['optional_intervals']
    observation['optional_intervals'][slug] = {'raw_entry_sequence': row['entry_sequence'],
                                             'raw_thread': row['thread'],
                                             'closed': False}
    optional_armed_stage = slug
    for breakpoint in optional_breakpoints:
        breakpoint.thread = row['thread']
        breakpoint.enabled = True


def close_optional_math(row):
    global optional_armed_stage
    slug = row['function'][len('sirius_simd_'):-len('_packet')]
    interval = observation['optional_intervals'][slug]
    interval.update(closed=True, raw_return_sequence=row['return_sequence'])
    if optional_armed_stage == slug:
        for breakpoint in optional_breakpoints:
            breakpoint.enabled = False
        optional_armed_stage = None


class OptionalReturned(gdb.FinishBreakpoint):
    def __init__(self, row):
        super().__init__(gdb.newest_frame(), internal=True)
        self.row = row

    def stop(self):
        try:
            row = self.row
            row.update(returned=True, completed=True, return_sequence=event_number(),
                       return_thread=gdb.selected_thread().global_num,
                       return_mxcsr=reg('mxcsr'))
            save()
            assert row['return_thread'] == row['thread']
            assert row['return_mxcsr'] & 0xe040 == 0
            return False
        except Exception as error:
            observation['failures'].append('optional return: ' + repr(error))
            save()
            return True


class OptionalEntered(gdb.Breakpoint):
    def __init__(self, name):
        self.name = name
        address = int(gdb.parse_and_eval('(void *)&' + name))
        super().__init__('*' + hex(address), internal=True)
        self.enabled = False

    def stop(self):
        global optional_armed_stage
        for breakpoint in optional_breakpoints:
            breakpoint.enabled = False
        expected_stage = optional_armed_stage
        optional_armed_stage = None
        row = {'function': self.name, 'entry_sequence': event_number(),
               'pc': reg('rip'), 'library': gdb.solib_name(reg('rip')),
               'pid': gdb.selected_inferior().pid,
               'thread': gdb.selected_thread().global_num,
               'mxcsr': reg('mxcsr'), 'returned': False,
               'expected_stage': expected_stage}
        observation['optional_math_calls'].append(row)
        save()
        try:
            assert row['library'] and Path(row['library']).resolve(strict=True) == provider
            row['provider_maps'] = [line for line in Path(f"/proc/{row['pid']}/maps").read_text().splitlines()
                                    if str(provider) in line]
            bases = {int(line.split()[0].split('-')[0], 16) - int(line.split()[2], 16)
                     for line in row['provider_maps']}
            base = row['pc'] - optional_symbols[self.name]
            assert base in bases
            row['provider_base'] = base
            cache = base + optional_symbols['__cpu_features2']
            row['cpu_features2_address'] = cache
            row['cpu_features2_first12'] = bytes(gdb.selected_inferior().read_memory(cache, 12)).hex()
            frames = []
            caller = gdb.newest_frame().older()
            for _ in range(16):
                if caller is None:
                    break
                frames.append({'name': caller.name(), 'pc': int(caller.pc())})
                caller = caller.older()
            row['caller_frames'] = frames
            matches = [(slug, frame) for frame in frames for slug, ranges in raw_ranges.items()
                       if any(lo <= frame['pc'] - base < hi or lo <= frame['pc'] - base - 1 < hi
                              for lo, hi in ranges)]
            stages = {slug for slug, _ in matches}
            assert stages == {expected_stage}, ('raw caller stage', stages, expected_stage)
            row['stage'] = expected_stage
            assert row['thread'] == observation['optional_intervals'][expected_stage]['raw_thread']
            row['raw_caller'] = matches[0][1]
            row['scratch'] = reg('rdi')
            if '_tagged_' in self.name:
                row.update(tags=reg('rsi'), destination=reg('edx') & 0xffffffff,
                           a_index=reg('ecx') & 0xffffffff,
                           b_index=reg('r8d') & 0xffffffff, mask=reg('r9d') & 0xffffffff)
            else:
                row.update(destination=reg('esi') & 0xffffffff,
                           a_index=reg('edx') & 0xffffffff,
                           b_index=reg('ecx') & 0xffffffff, mask=reg('r8d') & 0xffffffff)
            save()
            assert row['scratch'] and 0 < row['mask'] <= 255
            assert row['mxcsr'] & 0xe040 == 0
            assert bytes.fromhex(row['cpu_features2_first12'])[8] & 4
            OptionalReturned(row)
            return False
        except Exception as error:
            observation['failures'].append('optional entry: ' + repr(error))
            save()
            return True
