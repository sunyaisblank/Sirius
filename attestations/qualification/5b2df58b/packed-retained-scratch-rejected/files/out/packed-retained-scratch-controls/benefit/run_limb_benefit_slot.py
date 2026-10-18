"""Set child limits before one original frozen timestamp slot."""
from pathlib import Path
import os
import resource
import sys

sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
resource.setrlimit(resource.RLIMIT_FSIZE,(64*1024*1024,64*1024*1024))
root=Path.cwd()
assert len(sys.argv)==3
executable=root/sys.argv[1]
assert executable in [root/('bin/linux-gcc/packed-retained-scratch-controls/original_'+mode+'_controls')
                      for mode in ('baseline','candidate')]
assert sys.argv[2]=='--gtest_filter=RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults'
os.execv(str(executable),[str(executable),sys.argv[2]])
