"""Create only the fresh private arithmetic output/resource directory."""
from pathlib import Path
import shutil
R = Path.cwd().resolve()
W = Path(__file__).resolve().parent
T = R / 'bin/linux-gcc/computed-carry-exit-original-prerequisites/arithmetic'
assert not T.exists()
(T / 'resources/tests/backend').mkdir(parents=True)
shutil.copyfile(W / 'portable_binary32_reference.bin', T / 'resources/tests/backend/portable_binary32_reference.bin')
