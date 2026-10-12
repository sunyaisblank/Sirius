"""Run the unchanged original integrated filter with fixed process limits."""
from pathlib import Path
import os
import resource
import sys

sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
resource.setrlimit(resource.RLIMIT_FSIZE,(64*1024*1024,64*1024*1024))
root=Path.cwd()
executable=root/'bin/linux-gcc/portable-execution-attribution/limb-screen/original_candidate_controls'
assert len(sys.argv)==1
os.execv(str(executable),[str(executable),
    '--gtest_filter=RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults'])
