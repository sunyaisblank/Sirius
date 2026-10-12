"""Run the frozen original comparator once after the complete owned sequence."""
from pathlib import Path
import resource
import runpy
import sys

sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
resource.setrlimit(resource.RLIMIT_FSIZE,(64*1024*1024,64*1024*1024))
benefit=Path(__file__).resolve().parent/'benefit'
sys.path.insert(0,str(benefit))
runpy.run_path(str(benefit/'compare_software.py'),run_name='__main__')
