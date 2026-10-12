"""Passive JIT object collection for the single short-input compiler probe."""
import gdb
import hashlib
import json
import os
from pathlib import Path
import struct

folder = Path(os.environ['SIRIUS_LIMB_JIT_CAPTURE'])
events = []
state = {'inferior':None,'exit_code':None,'error':None,'events':events,
         'generated_function_breakpoints':False,'inferior_calls_injected':False}
total = 0
def save():
    (folder/'gdb-capture.json').write_text(json.dumps(state,indent=2)+'\n')
def memory(address,size):
    return bytes(gdb.selected_inferior().read_memory(address,size))
class Main(gdb.Breakpoint):
    def stop(self):
        pid = gdb.selected_inferior().pid
        stat = Path(f'/proc/{pid}/stat').read_text()
        state['inferior']={'pid':pid,'start_ticks':stat[stat.rfind(')')+2:].split()[19],
                           'boot_id':Path('/proc/sys/kernel/random/boot_id').read_text().strip()}
        save()
        self.enabled = False
        return False
class Jit(gdb.Breakpoint):
    def stop(self):
        global total
        try:
            descriptor=int(gdb.parse_and_eval('(unsigned long)&__jit_debug_descriptor'))
            version,action,entry,first=struct.unpack('<IIQQ',memory(descriptor,24))
            assert version==1 and action in (1,2) and len(events)<64
            item={'action':action,'entry':entry,'first_entry':first}
            if entry:
                following,previous,address,size=struct.unpack('<QQQQ',memory(entry,32))
                item.update(object_address=address,object_bytes=size)
                if action==1:
                    assert 0<size<32*1024*1024
                    total+=size
                    assert total<128*1024*1024
                    raw=memory(address,size)
                    assert raw[:6]==b'\x7fELF\x02\x01'
                    name=f'jit-{len(events):02d}.elf'
                    path=folder/name
                    assert not path.exists()
                    path.write_bytes(raw)
                    item.update(file=name,sha256=hashlib.sha256(raw).hexdigest())
            events.append(item)
            save()
            return False
        except Exception as error:
            state['error']=repr(error)
            save()
            return True
def exited(event):
    state['exit_code']=getattr(event,'exit_code',None)
    save()
Main('main',internal=True)
Jit('__jit_debug_register_code',internal=True)
gdb.events.exited.connect(exited)
save()
