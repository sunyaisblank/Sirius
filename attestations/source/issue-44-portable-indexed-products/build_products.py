"""Build all actual products into the ignored preset stage under fixed limits."""
from pathlib import Path
import resource
import subprocess
import sys

sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
resource.setrlimit(resource.RLIMIT_FSIZE,(64*1024*1024,64*1024*1024))
root=Path.cwd()
stage=root/'bin/linux-gcc/portable-indexed-limb-accumulation/products'
assert not stage.exists()
stage.mkdir()
tools=root/'bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools'
subprocess.run(['python3','-B',str(root/'scripts/build-retained-kernels.py'),
    '--output',str(stage/'retained_kernels.h'),'--compiler','/opt/slang/bin/slangc',
    '--assembler',str(tools/'spirv-as'),'--disassembler',str(tools/'spirv-dis'),
    '--validator',str(tools/'spirv-val'),'--optimizer','/usr/bin/spirv-opt'],check=True)
assert (stage/'retained_kernels.h').exists()
