"""Apply the byte-reviewed patch; do not inherit the disproved AJ tag claim."""
from pathlib import Path
base=Path.cwd()/'out/native-cpu-predicated-refinement/apply_provider.py'
source=base.read_text().replace('AK adds only empty-work original terminal selections in Divide/Sqrt; exact inverse restores AJ and any live refinement retains original arithmetic','AL makes twelve tag invalidation/refresh writes reachable before their original returns; AK arithmetic is unchanged and the original stale-tag counterexamples pass')
exec(compile(source,str(base),'exec'))
