"""Use the preset's exact adapter flags with isolated diagnostic output paths."""
from pathlib import Path
import json
import resource
import shlex
import subprocess
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024 * 1024, 64 * 1024 * 1024))
root = Path.cwd()
work = Path(__file__).resolve().parent
stage = root / 'bin/linux-gcc/portable-execution-attribution/limb-screen'
commands = json.loads((root/'bin/linux-gcc/compile_commands.json').read_text())
rows = [entry for entry in commands if entry['file'] == str(root/'src/sirius/backend/vulkan/vulkan_device.cpp')]
assert len(rows) == 1
original = shlex.split(rows[0]['command'])
common = []
at = 0
while at < len(original):
    if original[at] in ('-o', '-c'):
        assert at + 1 < len(original)
        at += 2
    else:
        common.append(original[at])
        at += 1
assert '-Werror' in common and '-Wno-missing-field-initializers' in common
assert '-DSIRIUS_RELEASE_RESOURCE_LOCKED=1' in common
obj = stage/'vulkan-device-v2.o'
executable = stage/'short_input_jit_probe_v2'
assert not obj.exists() and not executable.exists()
recipe = work/'short-input-preset-recipe.json'
assert not recipe.exists()
recipe.write_text(json.dumps({'recorded_CMake_adapter_command':rows[0],
                             'unchanged_compiler_flags':common,
                             'only_output_and_dependency_paths_replaced':True},indent=2)+'\n')
subprocess.run([*common,'-c',rows[0]['file'],'-o',str(obj),'-MMD','-MF',str(stage/'vulkan-device-v2.d')],check=True)
subprocess.run([*common,str(work/'short_input_jit_probe.cpp'),str(obj),
                str(root/'bin/linux-gcc/src/sirius/libsirius_base.a'),'-lvulkan','-ldl','-pthread',
                '-Wl,-Map='+str(stage/'short-input-v2.map'),'-o',str(executable),
                '-MMD','-MF',str(stage/'short-input-v2.d')],check=True)
assert executable.exists()
