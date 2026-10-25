from pathlib import Path
base=Path.cwd()/'out/native-cpu-slang-bridge/build_provider.py';source=base.read_text();anchor='native-cpu-slang-bridge/provider';assert source.count(anchor)==1
exec(compile(source.replace(anchor,'native-cpu-local-identities/provider'),str(base),'exec'))
