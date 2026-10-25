from pathlib import Path
base=Path.cwd()/'out/native-cpu-normalization-tags/apply_provider.py'
source=base.read_text().replace('AJ caches the exact private-register normalization predicate with complete write invalidation and unchanged arithmetic/fallback','AN performs original AM complete-fact identity stores in ISPC before the unchanged C++ helper; all general arithmetic, writers and fallback remain')
exec(compile(source,str(base),'exec'))
