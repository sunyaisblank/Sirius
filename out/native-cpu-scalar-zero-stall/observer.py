"""One first scalar Endpoint input capture and one fixed-time passive stack snapshot."""
from pathlib import Path
import gdb
import hashlib
import json
import struct

W=Path.cwd()/'out/native-cpu-scalar-zero-stall'
F=W/'capture'
record={'native_calls':[],'failures':[],'stack_snapshot':False,'inferior_exit':None,
    'scope':'Diagnostic only: original failed C3 process prefix, first scalar Endpoint ABI/input and one15s all-thread stack snapshot, then deliberate termination. No successful test/qualification, resampling or frame claim.'}


def save():
    (F/'observation.json').write_text(json.dumps(record,indent=2)+'\n')


def reg(name):return int(gdb.parse_and_eval('$'+name))


class Returned(gdb.FinishBreakpoint):
    def __init__(self,row):
        super().__init__(gdb.newest_frame(),internal=True);self.row=row
    def stop(self):
        self.row.update(returned=True,handled=bool(reg('rax')&0xff),return_mxcsr=reg('mxcsr'))
        save();return False


class Entered(gdb.Breakpoint):
    def __init__(self):
        super().__init__('*'+hex(int(gdb.parse_and_eval('(void *)&sirius_native_endpoint_try_packet'))),internal=True)
    def stop(self):
        self.enabled=False
        try:
            inferior=gdb.selected_inferior();stack=struct.unpack('<5Q',inferior.read_memory(reg('rsp')+8,40).tobytes())
            row={'function':'sirius_native_endpoint_try_packet','pid':inferior.pid,
                'thread':gdb.selected_thread().global_num,'pc':reg('rip'),'library':gdb.solib_name(reg('rip')),
                'input':reg('rdi'),'input_bytes':reg('esi')&0xffffffff,'output':reg('rdx'),'output_bytes':reg('ecx')&0xffffffff,
                'first_x':reg('r8d')&0xffffffff,'y':reg('r9d')&0xffffffff,'mxcsr':reg('mxcsr'),'returned':False}
            row.update(zip(('z','grid_x','grid_y','grid_z','groups'),(v&0xffffffff for v in stack)))
            record['native_calls'].append(row);save()
            assert 4<row['input_bytes']<=16*1024*1024
            data=inferior.read_memory(row['input'],row['input_bytes']).tobytes()
            with (F/'first-scalar-endpoint-input.bin').open('xb') as stream:stream.write(data)
            row['input_identity']={'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
            Returned(row);save();return False
        except Exception as error:
            record['failures'].append(repr(error));save();return True


armed=False
breakpoints=[]


def loaded(value):
    global armed
    if not armed and value.new_objfile.filename and Path(value.new_objfile.filename).name=='libvulkan_lvp.so':
        breakpoints.append(Entered());armed=True


def exited(value):
    record['inferior_exit']=getattr(value,'exit_code',None);save()


def capture_stacks():
    try:
        text=gdb.execute('thread apply all bt 16',to_string=True)
        (F/'all-thread-stacks.txt').write_text(text)
        (F/'shared-libraries.txt').write_text(gdb.execute('info sharedlibrary',to_string=True))
        record['stack_snapshot']=True
        record['inferior_threads']=len(gdb.selected_inferior().threads())
        save()
    except Exception as error:
        record['failures'].append(repr(error));save();raise


gdb.events.new_objfile.connect(loaded)
gdb.events.exited.connect(exited)
save()
