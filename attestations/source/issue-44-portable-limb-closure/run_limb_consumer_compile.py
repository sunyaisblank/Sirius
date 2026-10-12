"""Execute one unchanged preset-derived C++ recipe with bounded child output."""
from pathlib import Path
import json
import os
import resource
import sys

sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
resource.setrlimit(resource.RLIMIT_FSIZE,(64*1024*1024,64*1024*1024))
work=Path(__file__).resolve().parent
assert len(sys.argv)==2
action=json.loads((work/'actions.json').read_text())[sys.argv[1]]
command=action['compiler_argv']
assert command[0]=='/usr/bin/g++-14'
at=command.index('-o')
assert not Path(command[at+1]).exists()
os.execv(command[0],command)
