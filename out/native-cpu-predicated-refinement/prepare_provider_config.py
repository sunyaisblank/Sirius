"""Reuse current build/config closure; no raw/scalar/host compiler action."""
from pathlib import Path
base=Path.cwd()/'out/native-cpu-normalization-tags/prepare_provider_config.py'
source=base.read_text().replace('native-cpu-direct-compound/provider/actions.json','native-cpu-normalization-tags/provider/actions.json')
exec(compile(source,str(base),'exec'))
