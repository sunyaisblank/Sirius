"""Passively capture first actual rows/pure/field/layer calls and returns."""
from pathlib import Path
import gdb
import json
import struct

work = Path.cwd() / 'out/native-cpu-field-cohort'
action = json.loads((work / 'actions.json').read_text())['control-8']
folder = work / 'control-8'
provider = Path(action['runtime_provider']).resolve(strict=True)
observation = {'calls': [], 'preflight_calls': [], 'field_calls': [], 'layer_calls': [],
               'capability_calls': [], 'failures': [], 'inferior_exit': None}
breakpoints = []
armed = False
sequence = 0


def save():
    (folder / 'native-observation.json').write_text(json.dumps(observation, indent=2) + '\n')


def reg(name):
    return int(gdb.parse_and_eval('$' + name))


def event():
    global sequence
    sequence += 1
    return sequence


class Returned(gdb.FinishBreakpoint):
    def __init__(self, row):
        super().__init__(gdb.newest_frame(), internal=True)
        self.row = row

    def stop(self):
        row = self.row
        try:
            row.update(returned=True, return_sequence=event(), return_thread=gdb.selected_thread().global_num,
                       return_mxcsr=reg('mxcsr'), handled=bool(reg('rax') & 0xff))
            if row['category'] == 'layer_calls':
                row['status_after'] = struct.unpack('<I',gdb.selected_inferior().read_memory(row['status'],4).tobytes())[0]
            save()
            assert row['handled'], 'first observed entry refused: ' + row['function']
            assert row['return_thread'] == row['thread']
            if row['category'] in ('preflight_calls','capability_calls'):
                assert row['mxcsr'] == row['return_mxcsr']
            else:
                assert row['return_mxcsr'] & 0xe040 == 0
            return False
        except Exception as error:
            observation['failures'].append(row['function'] + ': ' + repr(error)); save()
            return True


class Entered(gdb.Breakpoint):
    def __init__(self, name, category):
        self.name, self.category = name, category
        super().__init__('*' + hex(int(gdb.parse_and_eval('(void *)&' + name))),internal=True)

    def stop(self):
        self.enabled = False
        try:
            inferior = gdb.selected_inferior()
            row = {'function':self.name,'category':self.category,'entry_sequence':event(),
                'pc':reg('rip'),'library':gdb.solib_name(reg('rip')),'pid':inferior.pid,
                'thread':gdb.selected_thread().global_num,'mxcsr':reg('mxcsr'),'returned':False}
            if self.category != 'capability_calls':
                row.update(input=reg('rdi'),output=reg('rdx'))
                row['capacity'] = struct.unpack('<I',inferior.read_memory(row['input'],4).tobytes())[0]
                if self.category == 'layer_calls':
                    stack = struct.unpack('<9Q',inferior.read_memory(reg('rsp')+8,72).tobytes())
                    row.update(input_words=reg('esi')&0xffffffff,output_words=reg('ecx')&0xffffffff,
                               program_base=reg('r8d')&0xffffffff,count=reg('r9d')&0xffffffff)
                    row.update(zip(('registers','layer','start','end','row','row_base'),(v&0xffffffff for v in stack[:6])))
                    row.update(zip(('scratch','coefficients','status'),stack[6:]))
                    row['status_before'] = struct.unpack('<I',inferior.read_memory(row['status'],4).tobytes())[0]
                    callers=[]; frame=gdb.newest_frame().older()
                    while frame and len(callers)<8:
                        callers.append({'function':frame.name(),'pc':int(frame.pc())}); frame=frame.older()
                    row['callers']=callers
                else:
                    stack=struct.unpack('<5Q',inferior.read_memory(reg('rsp')+8,40).tobytes())
                    row.update(input_bytes=reg('esi')&0xffffffff,output_bytes=reg('ecx')&0xffffffff,
                               first_x=reg('r8d')&0xffffffff,y=reg('r9d')&0xffffffff)
                    row.update(zip(('z','grid_x','grid_y','grid_z','groups'),(v&0xffffffff for v in stack)))
            observation[self.category].append(row); save()
            assert row['library'] and Path(row['library']).resolve(strict=True)==provider
            row['provider_maps']=[line for line in Path(f'/proc/{inferior.pid}/maps').read_text().splitlines() if str(provider) in line]
            assert row['provider_maps']
            if self.category not in ('preflight_calls','capability_calls'):
                assert row['mxcsr'] & 0xe040 == 0
            Returned(row); save()
            return False
        except Exception as error:
            observation['failures'].append(self.name + ': ' + repr(error)); save()
            return True


def loaded(value):
    global armed
    if armed or not value.new_objfile.filename or Path(value.new_objfile.filename).resolve()!=provider:
        return
    breakpoints.append(Entered('sirius_native_simd_supported','capability_calls'))
    for slug in ('transport','endpoint','dopri_phase'):
        for suffix,category in (('_rows_try_packet','calls'),('_admit_packet','preflight_calls'),('_field_try_packet','field_calls')):
            breakpoints.append(Entered('sirius_native_'+slug+suffix,category))
        breakpoints.append(Entered('sirius_simd_'+slug+'_field_layer','layer_calls'))
    armed=True


def exited(value):
    observation['inferior_exit']=getattr(value,'exit_code',None);save()


gdb.events.new_objfile.connect(loaded)
gdb.events.exited.connect(exited)
save()
