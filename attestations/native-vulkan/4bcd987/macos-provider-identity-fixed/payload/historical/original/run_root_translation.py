from pathlib import Path
import importlib.util,hashlib,runpy,sys
sys.dont_write_bytecode=True
work=Path(__file__).resolve().parent
assert not sys.flags.optimize
assert hashlib.sha256((work/'input-freeze.json').read_bytes()).hexdigest()=='71ca1547b3fb616ac99b9e8871a7e2410f25443131a6652b8b1f1c3464f526c9'
spec=importlib.util.spec_from_file_location('reader',work/'check_translation.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.frozen_inputs()
runpy.run_path(str(work/'run_owned.py'),run_name='__main__')
