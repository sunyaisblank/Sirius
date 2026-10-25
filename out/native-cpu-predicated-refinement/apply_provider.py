"""Reuse the exact three-file patch/source identity application checks."""
from pathlib import Path
base=Path.cwd()/'out/native-cpu-normalization-tags/apply_provider.py'
source=base.read_text().replace('AJ caches the exact private-register normalization predicate with complete write invalidation and unchanged arithmetic/fallback','AK adds only empty-work original terminal selections in Divide/Sqrt; exact inverse restores AJ and any live refinement retains original arithmetic')
exec(compile(source,str(base),'exec'))
