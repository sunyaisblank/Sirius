"""Reuse the original sole-C8 planner with current host/provider/component pins."""
from pathlib import Path

basis = Path.cwd() / 'out/native-cpu-zero-term-cohort/prepare.py'
source = basis.read_text().replace('out/native-cpu-zero-term-cohort', 'out/native-cpu-field-cohort')
source = source.replace('out/native-cpu-zero-term/', 'out/native-cpu-field-simd/')
source = source.replace('ispc-zero-compile-result.json', 'ispc-compile-result.json')
source = source.replace('nativeCPU providerS', 'nativeCPU providerT')
source = source.replace('unchanged Q Transport scheduling and new exact-zero3 objects; fresh Transport compatibility required',
    'field SIMD inside each original logical row task; fresh full38 Transport compatibility passed')
exec(compile(source, str(basis), 'exec'))
