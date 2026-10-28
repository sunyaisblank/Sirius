"""Read the existing AU owner, correcting only an inherited reader namespace.

The attempted workload and every frozen input remain unchanged. The original
reader wrapper translated its `current` evidence root to AT instead of AU;
this adapter corrects that exact generated literal before its checks execute.
"""
from pathlib import Path
import frame_authority

frame_authority.check_acceptance()
original_load = frame_authority.load


def load(name):
    source = original_load(name)
    if name == 'read_result.py':
        before = "current = root / 'out/native-cpu-dispatch-token-cohort'"
        after = "current = root / 'out/native-cpu-shared-value-cohort'"
        assert source.count(before) == 1
        source = source.replace(before, after, 1)
    return source


frame_authority.load = load
basis = Path(__file__).with_name('read_result.py')
exec(compile(basis.read_text(), str(basis), 'exec'),
     {'__file__': str(basis), '__name__': '__main__'})
