from pathlib import Path
base=Path.cwd()/'out/native-cpu-slang-bridge/prepare_provider_build.py'
exec(compile(base.read_text(),str(base),'exec'))
