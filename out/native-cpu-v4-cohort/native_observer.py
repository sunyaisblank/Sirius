"""Original C8 observation plus bounded first optional math per raw interval."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path.cwd() / 'out/native-cpu-v4-cohort'))
from rebind import load

source = load('native_observer.py')
extension = (Path.cwd() / 'out/native-cpu-v4-cohort/optional_observer.py').read_text()
anchor = 'def loaded(event):'
assert source.count(anchor) == 1
source = source.replace(anchor, extension + '\n\n' + anchor, 1)
anchor = "                row['completed'] = True  # Raw SIMD export returns void."
assert source.count(anchor) == 1
source = source.replace(anchor, anchor + '\n                close_optional_math(row)', 1)
anchor = '            Returned(row)\n            save()\n            return False'
assert source.count(anchor) == 3
before, after = source.rsplit(anchor, 1)
source = before + "            if self.category == 'simd_calls':\n                arm_optional_math(row)\n" + anchor + after
anchor = '    breakpoints.append(CapabilityEntered())'
assert source.count(anchor) == 1
source = source.replace(anchor, '''    for kind in ('dense', 'tagged'):
        for operation in ('add', 'multiply', 'subtract', 'divide', 'sqrt'):
            breakpoint = OptionalEntered('sirius_native_v4_' + kind + '_' + operation + '8')
            optional_breakpoints.append(breakpoint)
            breakpoints.append(breakpoint)
''' + anchor, 1)
exec(compile(source, str(Path.cwd() / 'out/native-cpu-local-identities-cohort/native_observer.py'), 'exec'))
