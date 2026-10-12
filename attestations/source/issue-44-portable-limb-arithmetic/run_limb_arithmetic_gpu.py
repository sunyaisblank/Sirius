"""Set fixed child limits before the ordinary synchronous Vulkan dispatch."""
from pathlib import Path
import os
import resource
import sys

sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
resource.setrlimit(resource.RLIMIT_FSIZE,(64*1024*1024,64*1024*1024))
root=Path.cwd()
work=Path(__file__).resolve().parent
build=root/'bin/linux-gcc/portable-execution-attribution/limb-screen'
assert len(sys.argv)==4
producer,mode,output=sys.argv[1:]
assert producer in ('baseline','candidate')
assert mode in ('integer','normal','normal_transport','normal_endpoint')
expected=work/'limb-arithmetic-capture-shape.bin' if producer=='baseline' else work/('gpu-limb-arithmetic-baseline-'+mode)/'actual.bin'
executable=build/'limb_arithmetic_gpu'
os.execv(str(executable),[str(executable),str(build/('arithmetic-'+producer+'-'+mode+'.spv')),
    str(work/'limb-arithmetic-input.bin'),str(expected),output,
    '--capture' if producer=='baseline' else '--dispatch'])
