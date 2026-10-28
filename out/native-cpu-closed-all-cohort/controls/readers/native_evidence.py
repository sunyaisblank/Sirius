from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from control_authority import ao
exec(compile(ao('readers/native_evidence.py'), '/home/astra/.project/Sirius/out/native-cpu-v4-cohort/controls/readers/native_evidence.py', 'exec'))
