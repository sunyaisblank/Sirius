"""Inspect the second actual pair already observed in the interval inventory."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
WORK=Path(__file__).resolve().parent
ps=(WORK/'profiler_trace_samples.ps1').read_text().replace('profiler-trace-samples','profiler-interval-second-pair').replace('4 vkCmdDispatch(3, 1, 1)','16 vkCmdDispatch(3, 1, 1)').replace('10 vkCmdDispatch(6, 1, 1)','22 vkCmdDispatch(6, 1, 1)')
b=(WORK/'run_profiler_trace_samples.py').read_text().replace('profiler-trace-samples','profiler-interval-second-pair').replace('profiler_trace_samples.ps1','profiler_interval_second_pair.ps1').replace('trace-samples-source-review.json','interval-second-pair-source-review.json').replace('prepare_trace_samples.py','prepare_interval_pair2.py').replace('native-capture-early-console','native-capture-interval-console').replace('original-legacy-te.rgp','original-legacy-interval.rgp').replace('4 vkCmdDispatch(3, 1, 1)','16 vkCmdDispatch(3, 1, 1)').replace('10 vkCmdDispatch(6, 1, 1)','22 vkCmdDispatch(6, 1, 1)')
old="WORK/'profiler-trace-owned-click/transport-automation.json'"
assert b.count(old)==1
b=b.replace(old,old+",WORK/'profiler-interval-inventory-bridge.json',WORK/'profiler-interval-inventory/automation.json'")
for name in ('profiler_interval_second_pair.ps1','run_profiler_interval_second_pair.py'):assert not (WORK/name).exists()
(WORK/'profiler_interval_second_pair.ps1').write_text(ps)
(WORK/'run_profiler_interval_second_pair.py').write_text(b)
