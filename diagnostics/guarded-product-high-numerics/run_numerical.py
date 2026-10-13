"""Set fixed limits and execute the unchanged raw-word GPU comparison path."""
from pathlib import Path
import os
import resource
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024**2, 64 * 1024**2))
root = Path.cwd().resolve()
work = Path(__file__).resolve().parent
stage = root / 'bin/linux-gcc/guarded-product-high-numerics'
family, mode, output = sys.argv[1:]
assert family in ['residual', 'limb']
assert mode in (['core', 'rte32'] if family == 'residual' else
               ['integer', 'normal', 'normal_transport', 'normal_endpoint'])
path = Path(output).resolve()
assert path.is_relative_to(work) and not path.exists()
executable = stage / ('residual_gpu' if family == 'residual' else 'limb_arithmetic_gpu')
fixture = work / ('residual-input.bin' if family == 'residual' else 'limb-arithmetic-input.bin')
expected = work / 'residual-expected.bin' if family == 'residual' else stage / 'baseline-outputs' / f'{mode}.bin'
module = stage / f'{family}-candidate-{mode}.spv'
os.execv(str(executable), [str(executable), str(module), str(fixture), str(expected), str(path), '--dispatch'])
