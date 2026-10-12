"""Use unchanged preset flags and the already fresh exact Vulkan adapter."""
from pathlib import Path
import json
import resource
import subprocess
import sys

sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
resource.setrlimit(resource.RLIMIT_FSIZE,(64*1024*1024,64*1024*1024))
root=Path.cwd()
work=Path(__file__).resolve().parent
build=root/'bin/linux-gcc/portable-execution-attribution/limb-screen'
recipe=json.loads((work/'short-input-preset-recipe.json').read_text())
common=recipe['unchanged_compiler_flags']
assert '-Werror' in common and '-Wno-missing-field-initializers' in common
output=build/'limb_arithmetic_gpu'
assert not output.exists()
subprocess.run([*common,str(work/'limb_arithmetic_gpu.cpp'),str(build/'vulkan-device-v2.o'),
    str(root/'bin/linux-gcc/src/sirius/libsirius_base.a'),'-lvulkan','-ldl','-pthread',
    '-Wl,-Map='+str(build/'limb-arithmetic-gpu.map'),'-o',str(output),
    '-MMD','-MF',str(build/'limb-arithmetic-gpu.d')],check=True)
