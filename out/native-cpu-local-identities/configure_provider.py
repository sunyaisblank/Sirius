from pathlib import Path
base=Path.cwd()/'out/native-cpu-slang-bridge/configure_provider.py'
exec(compile(base.read_text(),str(base),'exec'))
