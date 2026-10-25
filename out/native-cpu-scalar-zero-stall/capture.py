"""One bounded15s failure-stage observation; deliberately kill after the snapshot."""
from pathlib import Path
import json
import signal
import subprocess
import sys
import time

R=Path.cwd();W=Path(__file__).resolve().parent;F=W/'capture'
assert sys.argv[1:]==['capture',str(F)]
action=json.loads((W/'actions.json').read_text())['capture']
commands=F/'commands.gdb'
commands.write_text('set pagination off\nset confirm off\nset debuginfod enabled off\nset auto-load off\n'
    'set startup-with-shell off\nset disable-randomization off\nset print thread-events off\n'
    'source '+str(W/'observer.py')+'\nrun\npython capture_stacks()\nkill\nquit\n')
argv=['/usr/bin/gdb','-nx','-nh','-q','-batch','-x',str(commands),'--args',
      action['executable'],'--gtest_filter='+action['case'],'--gtest_color=no']
child=subprocess.Popen(argv,cwd=R)
started=time.monotonic();requested=False
while child.poll() is None:
    if not requested and time.monotonic()-started>=15:
        # This is our unreaped child; its PID cannot be reused.
        (F/'interrupt.json').write_text(json.dumps({'debugger_pid':child.pid,'elapsed_seconds':time.monotonic()-started,
            'signal_target':'own unreaped GDB child','signals_sent':1})+'\n')
        child.send_signal(signal.SIGINT);requested=True
    time.sleep(.05)
assert child.returncode==0,child.returncode
record=json.loads((F/'observation.json').read_text())
assert requested and record['stack_snapshot'] and not record['failures'],record
print(json.dumps({'diagnostic_stack_captured':True,'original_C3_test_accepted':False,
    'first_scalar_calls_captured':len(record['native_calls']),'scope':record['scope']}))
