"""One passive observation of the fixed positive case; loaded by positive.gdb."""
from pathlib import Path
import gdb
import hashlib
import importlib.util
import json
import sys
sys.dont_write_bytecode=True

root=Path.cwd()
work=root/'out/native-cpu-exact-product-residual'
executable=Path(gdb.current_progspace().filename).resolve()
target='sirius_native_v4_dense_multiply8'
convert_offset,post_offset=0x1e6a,0x1e70


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())


# Bind this one actual linked instruction to the already inspected object.
assert executable==root/'bin/linux-gcc/native-cpu-exact-product-residual/positive-check'
helper=root/'out/native-cpu-accumulate-compatibility/elf_proof.py'
assert pin(helper)['sha256']=='278766535b0ddcb9ad2ab97b68805cd851dda31036567a35f8bf5a07afedb0f7'
spec=importlib.util.spec_from_file_location('positive_existing_elf',helper)
elf_module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(elf_module)
linked=elf_module.Elf(executable)
symbol,=[row for row in linked.symbols if row['name']==target and row['section'] and row['size']]
code=linked.function_bytes(symbol)
assert code[convert_offset-6:post_offset+6]==bytes.fromhex(
    '62a1d5405cee62a1fd485aed62b15528dbed') # subpd; cvtpd2ps ymm21; andd.
assert code[0x1e86:0x1e88]==bytes.fromhex('7465') # All-nonzero jumps to+0x1eed.
compiled=json.loads((work/'compile-result.json').read_text())
optional=root/'bin/linux-gcc/native-cpu-exact-product-residual/optional.o'
assert compiled['exit']==0 and pin(optional)==compiled['optional_object']
obj=elf_module.Elf(optional)
original,=[row for row in obj.symbols if row['name']==target and row['section'] and row['size']]
a,b=bytearray(obj.function_bytes(original)),bytearray(code)
assert len(a)==len(b)
for section,location,kind,_addend,_target in obj.relocations:
    offset=location-original['value']
    if section==original['section'] and 0<=offset<original['size']:
        assert kind in (2,4,42) and offset+4<=original['size']
        a[offset:offset+4]=b[offset:offset+4]=b'\0'*4
assert a==b
record=dict(executable=dict(path=str(executable),**pin(executable)),optional_object=pin(optional),
    function=target,function_address=symbol['value'],function_bytes=symbol['size'],
    original_object_function_opcodes_match=True,convert_offset=convert_offset,
    after_convert_offset=post_offset,residual_register='ymm21',
    main_entry=None,residual=None,exit_code=None,coverage_passed=False)


def words(name):
    return [int(gdb.parse_and_eval('$'+name+'.v8_int32['+str(i)+']'))&0xffffffff for i in range(8)]


class Residual(gdb.Breakpoint):
    def __init__(self,base):
        super().__init__('*'+hex(base+post_offset),internal=True)
        self.base=base
    def stop(self):
        self.enabled=False
        row=dict(pc=int(gdb.parse_and_eval('$pc')),function_base=self.base,
            thread=gdb.selected_thread().global_num,mxcsr=int(gdb.parse_and_eval('$mxcsr')),
            operand_a=words('ymm3'),operand_b=words('ymm0'),high=words('ymm24'),residual=words('ymm21'))
        row['source_guard_true']=all(24<=(a&0x7fffffff)>>23 and 24<=(b&0x7fffffff)>>23 and
            (a&0x7fffffff)<=0x5d800000 and (b&0x7fffffff)<=0x5d800000 and
            ((a&0x7fffffff)>>23)+((b&0x7fffffff)>>23)>=174
            for a,b in zip(row['operand_a'],row['operand_b']))
        record['residual']=row
        return False


class Main(gdb.Breakpoint):
    def __init__(self):
        super().__init__('main',internal=True)
    def stop(self):
        self.enabled=False
        base=int(gdb.parse_and_eval('(void*)&'+target))
        assert bytes(gdb.selected_inferior().read_memory(base+convert_offset,6)).hex()=='62a1fd485aed'
        cache=int(gdb.parse_and_eval('(void*)&__cpu_features2'))
        data=bytes(gdb.selected_inferior().read_memory(cache,12))
        record['main_entry']=dict(pc=int(gdb.parse_and_eval('$pc')),function_base=base,
            load_base=base-symbol['value'],cpu_features2_address=cache,
            cpu_features2_first12=data.hex(),actual_v4_cache_bit=bool(data[8]&4))
        Residual(base)
        return False


Main()
def exited(event):
    record['exit_code']=getattr(event,'exit_code',None)
gdb.events.exited.connect(exited)


def finish():
    row,entry=record['residual'],record['main_entry']
    record['coverage_passed']=bool(entry and row and record['exit_code']==0 and
        entry['actual_v4_cache_bit'] and row['source_guard_true'] and
        row['operand_a']==row['operand_b']==[0x3f800001]*8 and
        row['high']==[0x3f800002]*8 and row['residual']==[0x28800000]*8 and
        row['pc']==entry['function_base']+post_offset and row['mxcsr']&0xe040==0)
    with (work/'positive-run/targets.json').open('x') as stream:
        stream.write(json.dumps(record,indent=2)+'\n')
    if not record['coverage_passed']:
        gdb.execute('quit 1')
