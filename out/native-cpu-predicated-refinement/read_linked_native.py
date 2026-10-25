"""Reuse all actual AJ object/code/module/relocated-call joins."""
from pathlib import Path
base=Path.cwd()/'out/native-cpu-normalization-tags/read_linked_native.py'
exec(compile(base.read_text(),str(base),'exec'))
