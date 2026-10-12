"""Compile the real portable seed paths with the production SSA preparation."""
from pathlib import Path
import subprocess
import sys

root = Path.cwd()
work = Path(__file__).resolve().parent
build = root / 'bin/linux-gcc/portable-denominator-preparation'
assert len(sys.argv) == 2 and sys.argv[1] in ('baseline', 'candidate')
producer = sys.argv[1]
include = build / 'baseline' if producer == 'baseline' else root / 'src/sirius/kernels'
definitions = ['-DSIRIUS_RETAINED_PORTABLE=1']
if producer == 'baseline':
    definitions.append('-DSIRIUS_PREPARED_BASELINE=1')
raw = build / (producer + '-seed.compiler.spv')
module = build / (producer + '-seed.spv')
assembly = build / (producer + '-seed.spvasm')
tool = root / 'bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools'
subprocess.run(['/opt/slang/bin/slangc', str(work / 'seed_probe.slang'), *definitions,
                '-I', str(include), '-O0', '-target', 'spirv', '-profile', 'spirv_1_5',
                '-entry', 'ComputeMain', '-stage', 'compute', '-o', str(raw)], check=True)
subprocess.run(['/usr/bin/spirv-opt', '--target-env=vulkan1.2', '--ssa-rewrite',
                '--eliminate-dead-code-aggressive', '--preserve-bindings', '--preserve-interface',
                str(raw), '-o', str(module)], check=True)
subprocess.run([str(tool / 'spirv-val'), '--target-env', 'vulkan1.2', str(module)], check=True)
subprocess.run([str(tool / 'spirv-dis'), str(module), '-o', str(assembly)], check=True)
