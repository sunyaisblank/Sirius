"""Fresh bounded device adapter and diagnostic; no canonical product build."""
from pathlib import Path
import resource
import subprocess
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024 * 1024, 64 * 1024 * 1024))
root = Path.cwd()
work = Path(__file__).resolve().parent
stage = root / 'bin/linux-gcc/portable-execution-attribution/limb-screen'
obj = stage / 'vulkan-device.o'
executable = stage / 'short_input_jit_probe'
assert not obj.exists() and not executable.exists()
common = ['/usr/bin/g++-14', '-std=c++2c', '-O3', '-DNDEBUG', '-Wall', '-Wextra',
          '-Wpedantic', '-Werror', '-I'+str(root/'src'), '-DSIRIUS_CONTRACT_MODE=2',
          '-DSIRIUS_HAS_RETAINED_COMPUTE=1']
subprocess.run([*common, '-c', str(root/'src/sirius/backend/vulkan/vulkan_device.cpp'),
                '-o', str(obj), '-MMD', '-MF', str(stage/'vulkan-device.d')], check=True)
subprocess.run([*common, str(work/'short_input_jit_probe.cpp'), str(obj),
                str(root/'bin/linux-gcc/src/sirius/libsirius_base.a'), '-lvulkan', '-ldl', '-pthread',
                '-Wl,-Map='+str(stage/'short-input.map'), '-o', str(executable),
                '-MMD', '-MF', str(stage/'short-input.d')], check=True)
assert executable.exists()
