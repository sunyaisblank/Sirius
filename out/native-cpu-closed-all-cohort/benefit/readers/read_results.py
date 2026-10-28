"""Retain every original numeric/provenance gate with current AX composition scope."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from authority import load, an, source as pinned_source, reader_scope

load('readers/read_results.py')
an('readers/read_results.py')
basis = Path.cwd() / 'out/native-cpu-active-prefix/benefit/readers/read_results.py'
source = pinned_source(basis, '1371748434db9de22af3555de6e2b27b25e9bd25f7928d1937f8f0202a7ace92')
assert 'out/native-cpu-active-prefix' in source
source = source.replace('out/native-cpu-active-prefix', 'out/native-cpu-closed-all-cohort')
anchor = "exec(compile(source, str(basis), 'exec'))"
assert source.count(anchor) == 1
source = source.replace(anchor, 'source = reader_scope(source)\n' + anchor, 1)
exec(compile(source, str(basis), 'exec'))
