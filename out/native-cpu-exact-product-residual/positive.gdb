set pagination off
set confirm off
set debuginfod enabled off
set auto-load off
set disable-randomization off
set startup-with-shell off
python
from pathlib import Path
path=Path.cwd()/'out/native-cpu-exact-product-residual/positive.py'
exec(compile(path.read_text(),str(path),'exec'))
end
run
python
finish()
end
