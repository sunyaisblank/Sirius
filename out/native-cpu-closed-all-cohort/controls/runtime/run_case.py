"""Retain original ABBA and preceding-control prerequisites for each case."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from control_authority import ao, load

ao('runtime/run_case.py')
exec(compile(load('runtime/run_case.py'), 'out/native-cpu-local-identities-cohort/controls/runtime/run_case.py', 'exec'))
