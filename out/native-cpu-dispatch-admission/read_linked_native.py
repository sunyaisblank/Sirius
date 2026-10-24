"""Preserve exact ELF associations; follow the actual shared guard call chain."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-slang-bridge/read_linked_native.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '9a94c2b45471807ad74b28060edc60e19f33d42e7b9f94ba6126065d79c6f064'
source = basis.read_text()
changes = [
    ("batch = json.loads((W / 'batch-generation.json').read_text())",
     "batch = json.loads((R / 'out/native-cpu-slang-bridge/batch-generation.json').read_text())"),
    ("disassembly = subprocess.check_output(['objdump', '-d', '--disassemble=cs_exec_fn', str(library)], text=True)",
     "callback = subprocess.check_output(['objdump', '-d', '--disassemble=cs_exec_fn', str(library)], text=True)\n"
     "assert '<lp_native_cpu_packet_mode>' in callback\n"
     "packet_mode = subprocess.check_output(['objdump', '-d', '--disassemble=lp_native_cpu_packet_mode', str(library)], text=True)\n"
     "for slug in ('transport', 'endpoint', 'dopri_phase'):\n"
     "    assert '<sirius_native_' + slug + '_admit_packet>' in packet_mode\n"
     "disassembly = callback + packet_mode"),
    ("'callback_direct_calls_all3': True", "'callback_via_shared_packet_mode_reaches_all3': True"),
    ('Callback includes direct native calls and ordinary fallback under existing IEEE32 guard.',
     'Callback reaches all three native executors through the shared packet-mode helper under existing IEEE32 guard, with ordinary fallback. The helper also directly reaches all three pure admission predicates.'),
]
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
