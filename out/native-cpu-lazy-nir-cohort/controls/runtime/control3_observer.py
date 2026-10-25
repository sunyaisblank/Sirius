"""Observe actual finalizer/fallback routes during the one original V C3 run.

No extra dispatch, fixed-time sample, inferior call or state write. Original
numerical assertions and300s/RSS/reserve limits remain authoritative.
"""
from pathlib import Path
import gdb
import json
import time

W=Path.cwd()/'out/native-cpu-lazy-nir-cohort/controls/runtime'
F=W/'control-3'
action=json.loads((W/'actions.json').read_text())['control-3']
provider=Path(action['runtime_provider']).resolve(strict=True)
layout=json.loads((W/'nir-layout.json').read_text())
record={'compiler_calls':[],'native_calls':[],'failures':[],'inferior_exit':None,
    'scope':'Passive first64 compiler entries/returns and first Endpoint rows entry/return in original C3. Metadata offsets come from exact current finalizer code. No separate profile, sample, dispatch, retry, timing benefit or stronger numerical gate.'}
breakpoints=[];armed=False;sequence=0


def save():(F/'route-observation.json').write_text(json.dumps(record,indent=2)+'\n')


def event():
    global sequence
    sequence+=1;return sequence


def reg(name):return int(gdb.parse_and_eval('$'+name))


class Returned(gdb.FinishBreakpoint):
    def __init__(self,row):
        super().__init__(gdb.newest_frame(),internal=True);self.row=row
    def stop(self):
        row=self.row
        row.update(returned=True,return_sequence=event(),return_thread=gdb.selected_thread().global_num,
                   return_monotonic=time.monotonic())
        if row['category']=='native_calls':row['handled']=bool(reg('rax')&0xff)
        elif row['kind']=='finalize':row['returned_error_pointer']=reg('rax')
        save();return False


class Entered(gdb.Breakpoint):
    def __init__(self,name,kind):
        self.name,self.kind=name,kind
        super().__init__('*'+hex(int(gdb.parse_and_eval('(void *)&'+name))),internal=True)
    def stop(self):
        if self.kind=='native':self.enabled=False
        elif len(record['compiler_calls'])>=64:self.enabled=False;return False
        try:
            inferior=gdb.selected_inferior();pc=reg('rip');library=gdb.solib_name(pc)
            row={'function':self.name,'kind':self.kind,'pid':inferior.pid,
                'thread':gdb.selected_thread().global_num,'pc':pc,'library':library,
                'entry_sequence':event(),'entry_monotonic':time.monotonic(),'returned':False,
                'category':'native_calls' if self.kind=='native' else 'compiler_calls'}
            caller=gdb.newest_frame().older()
            row['caller_function']=caller.name() if caller else None
            if self.kind!='native':
                nir=reg('rsi' if self.kind=='finalize' else 'rdi');row['nir']=nir
                row['shader_stage']=inferior.read_memory(nir+layout['shader_stage_offset'],1).tobytes()[0]
                row['native_stage']=inferior.read_memory(nir+layout['native_stage_offset'],1).tobytes()[0]
            record[row['category']].append(row);save()
            assert library and Path(library).resolve(strict=True)==provider
            Returned(row);return False
        except Exception as error:
            record['failures'].append(self.name+': '+repr(error));save();return True


def loaded(value):
    global armed
    if armed or not value.new_objfile.filename or Path(value.new_objfile.filename).resolve()!=provider:return
    for name,kind in (('llvmpipe_finalize_nir','finalize'),('lp_build_opt_nir','optimize'),
                      ('sirius_native_endpoint_rows_try_packet','native')):
        breakpoints.append(Entered(name,kind))
    armed=True


def exited(value):record['inferior_exit']=getattr(value,'exit_code',None);save()


gdb.events.new_objfile.connect(loaded)
gdb.events.exited.connect(exited)
save()
