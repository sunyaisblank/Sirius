"""Reuse the exercised sealed build planner; compile only the changed helper."""
from pathlib import Path
import hashlib
R=Path.cwd();base=R/'out/native-cpu-normalization-tags/prepare_cpp.py';source=base.read_text()
source=source.replace('native-cpu-normalization-tags','native-cpu-predicated-refinement')
source=source.replace('native-cpu-direct-compound-cohort','native-cpu-normalization-tags-cohort')
source=source.replace('native-cpu-direct-compound-frame-timeout','native-cpu-normalization-tags-frame-timeout')
source=source.replace("for slug in ('retained8', 'driver'):","for slug in ('retained8',):")
source=source.replace("str(B / 'driver.o')","str(R / 'bin/linux-gcc/native-cpu-normalization-tags/probe/driver.o')")
anchor="    compiler = Path(flags[0]).resolve(strict=True)"
extra="""    reused_driver=R/'bin/linux-gcc/native-cpu-normalization-tags/probe/driver.o'
    previous=doc(R/'out/native-cpu-normalization-tags/cpp-build-result.json')
    driver,=[r['object'] for r in previous['compilations'] if r['name']=='driver']
    assert pin(reused_driver)=={k:driver[k] for k in ('bytes','sha256')}
    for path in (reused_driver,R/'out/native-cpu-normalization-tags/driver.cpp',R/'out/native-cpu-normalization-tags/cpp-build-result.json'):
        inputs[str(path.relative_to(R))]=pin(path)
"""
assert source.count(anchor)==1;source=source.replace(anchor,extra+anchor)
source=source.replace("built['new_CPP_compilations'] == 2","built['new_CPP_compilations'] == 1")
start=source.index("                'scope': 'One new strict-FP helper")
end=source.index("}\n",start)
source=source[:start]+"                'scope': 'One changed strict-FP C++ helper restores six empty-work terminal returns in Divide/Sqrt, with exact inverse to AJ. Reuse its exact original 49-case driver object, original I adapter and all raw/scalar objects. Any nonempty work mask retains original seed, three ordered refinements, final residual/radius/refusal. Complete whole-bank/tag/mask/alias checks remain. No FP/ISA/count/cache/governor or scientific/release claim.'"+source[end:]
exec(compile(source,str(base),'exec'))
