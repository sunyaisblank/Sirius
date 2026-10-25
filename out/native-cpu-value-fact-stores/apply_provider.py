from pathlib import Path
base=Path.cwd()/'out/native-cpu-normalization-tags/apply_provider.py'
source=base.read_text().replace('AJ caches the exact private-register normalization predicate with complete write invalidation and unchanged arithmetic/fallback','AM caches complete exact value facts in the existing private word, with reachable invalidation/refresh and original identity-store precedence; all general arithmetic/fallback remains')
exec(compile(source,str(base),'exec'))
