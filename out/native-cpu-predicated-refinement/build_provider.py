from pathlib import Path
base=Path.cwd()/'out/native-cpu-slang-bridge/build_provider.py'
source=base.read_text();anchor='native-cpu-slang-bridge/provider';assert source.count(anchor)==1
source=source.replace(anchor,'native-cpu-predicated-refinement/provider')
exec(compile(source,str(base),'exec'))
