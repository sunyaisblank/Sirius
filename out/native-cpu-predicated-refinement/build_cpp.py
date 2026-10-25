"""Reuse original helper builder with one prepared compilation and old driver."""
from pathlib import Path
base=Path.cwd()/'out/native-cpu-normalization-tags/build_cpp.py'
source=base.read_text().replace("'new_CPP_compilations': 2, 'new_links': 1","'new_CPP_compilations': 1, 'new_links': 1")
exec(compile(source,str(base),'exec'))
