"""Keep the exact sole A1/B1/B2/A2 authority after current AR C8/static gates."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from authority import load
from correctness import check

check()
exec(compile(load('readers/run_abba.py'), 'out/native-cpu-v4-cohort/benefit/readers/run_abba.py', 'exec'))
