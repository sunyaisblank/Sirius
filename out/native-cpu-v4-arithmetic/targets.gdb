set pagination off
set confirm off
set startup-with-shell off
python
import gdb,json
from pathlib import Path
record={'v4_entries':{},'v3_entries':{},'exit_code':None,
        'cpu_features2_first12':None}
class First(gdb.Breakpoint):
    def __init__(self,name,kind):
        super().__init__(name,internal=True)
        self.target_name=name
        self.kind=kind
    def stop(self):
        if record['cpu_features2_first12'] is None:
            address=int(gdb.parse_and_eval('(void*)&__cpu_features2'))
            record['cpu_features2_first12']=bytes(gdb.selected_inferior().read_memory(address,12)).hex()
        record[self.kind][self.target_name]={'pc':int(gdb.parse_and_eval('$pc')),
            'thread':gdb.selected_thread().global_num}
        self.enabled=False
        return False
for version in ('v3','v4'):
    for operation in ('add','multiply','divide','sqrt','subtract'):
        First('sirius_native_'+version+'_dense_'+operation+'8',version+'_entries')
def exited(event):
    record['exit_code']=getattr(event,'exit_code',None)
gdb.events.exited.connect(exited)
end
run
python
path=Path.cwd()/'out/native-cpu-v4-arithmetic/dense-run/targets.json'
with path.open('x') as stream:
    stream.write(json.dumps(record,indent=2)+'\n')
if record['exit_code']!=0:
    gdb.execute('quit 1')
end
