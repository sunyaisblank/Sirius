"""Select the two observed owned tabs to inspect saved trace event timing."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
WORK = Path(__file__).resolve().parent
def replace(s, old, new):
    assert s.count(old) == 1, (old, s.count(old))
    return s.replace(old, new)
def save(n, s):
    p = WORK / n
    assert not p.exists()
    p.write_text(s)
ps = (WORK / 'profiler_trace_events.ps1').read_text().replace(
    'profiler-trace-events', 'profiler-trace-timing')
start = ps.index('    $tabName=')
end = ps.index('    $record.navigation_completed=$true', start)
block = ps[start:end]
block = block.replace("::NameProperty,'EVENTS'", '::NameProperty,$tabLabel')
block = replace(block, "$record.navigation_requested='EVENTS SelectionItemPattern.Select on exact owned root'", "$record.navigation_requested=$tabLabel+' SelectionItemPattern.Select on exact owned root'")
block += "    $record.navigation_history+=@($tabLabel)\n"
ps = ps[:start] + "    $record.navigation_history=@()\n    foreach ($tabLabel in @('EVENTS','Event timing')) {\n" + block + "    }\n" + ps[end:]
ps = replace(ps, 'Owned sealed-profile EVENTS tab selection and GUI observation only;', 'Owned sealed-profile EVENTS and Event timing tab selection and GUI observation only;')
ps = replace(ps, 'Loaded profile readiness lost after EVENTS selection', 'Loaded profile readiness lost after Event timing selection')
save('profiler_trace_timing.ps1', ps)
b = (WORK / 'run_profiler_trace_events.py').read_text().replace(
    'profiler-trace-events', 'profiler-trace-timing').replace(
    'profiler_trace_events.ps1', 'profiler_trace_timing.ps1').replace(
    'trace-events-source-review.json', 'trace-timing-source-review.json').replace(
    'prepare_trace_events.py', 'prepare_trace_timing.py')
b = replace(b, "o['navigation_completed'] and", "o['navigation_history']==['EVENTS','Event timing'] and o['navigation_completed'] and")
b = replace(b, 'Owned saved-profile EVENTS tab selection and window observation only.', 'Owned saved-profile EVENTS and Event timing selection and window observation only.')
b = replace(b, "WORK/'profiler-trace-stable/automation.json'", "WORK/'profiler-trace-stable/automation.json',WORK/'profiler-trace-events-bridge.json',WORK/'profiler-trace-events/automation.json'")
save('run_profiler_trace_timing.py', b)
