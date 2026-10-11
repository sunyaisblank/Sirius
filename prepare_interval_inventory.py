"""Read the acquired interval event layout before imposing a source join."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
WORK = Path(__file__).resolve().parent
ps=(WORK/'profiler_trace_timing.ps1').read_text().replace('profiler-trace-timing','profiler-interval-inventory')
b=(WORK/'run_profiler_trace_timing.py').read_text().replace('profiler-trace-timing','profiler-interval-inventory').replace('profiler_trace_timing.ps1','profiler_interval_inventory.ps1').replace('trace-timing-source-review.json','interval-inventory-source-review.json').replace('prepare_trace_timing.py','prepare_interval_inventory.py').replace('native-capture-early-console','native-capture-interval-console').replace('original-legacy-te.rgp','original-legacy-interval.rgp')
for name in ('profiler_interval_inventory.ps1','run_profiler_interval_inventory.py'):assert not (WORK/name).exists()
(WORK/'profiler_interval_inventory.ps1').write_text(ps)
(WORK/'run_profiler_interval_inventory.py').write_text(b)
