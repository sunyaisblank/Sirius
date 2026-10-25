from pathlib import Path
base=Path.cwd()/'out/native-cpu-predicated-refinement/prepare_provider_config.py'
exec(compile(base.read_text(),str(base),'exec'))
