"""Preserve original passive worker capture; also observe first pure admissions."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-active-prefix/native_observer.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'c38d87c7cea09dd50b53750982987ab5bfdd46859b53d0d61598a9e7c66377fa'
source = basis.read_text()
changes = [
    ('out/native-cpu-active-prefix', 'out/native-cpu-dispatch-admission'),
    ("observation = {'calls': [], 'failures': [], 'inferior_exit': None}",
     "observation = {'calls': [], 'preflight_calls': [], 'failures': [], 'inferior_exit': None}"),
    ("            assert self.row['return_mxcsr'] & 0xe040 == 0, 'return IEEE32 execution state mismatch'",
     "            if self.row['is_preflight']:\n                assert self.row['return_mxcsr'] == self.row['mxcsr'], 'pure admission changed FP state'\n            else:\n                assert self.row['return_mxcsr'] & 0xe040 == 0, 'return IEEE32 execution state mismatch'"),
    ('    def __init__(self, name):\n        self.name = name',
     '    def __init__(self, name, is_preflight=False):\n        self.name = name\n        self.is_preflight = is_preflight'),
    ("'mxcsr': reg('mxcsr'), 'returned': False}",
     "'mxcsr': reg('mxcsr'), 'returned': False, 'is_preflight': self.is_preflight}"),
    ("            observation['calls'].append(row)",
     "            observation['preflight_calls' if self.is_preflight else 'calls'].append(row)"),
    ("            assert row['mxcsr'] & 0xe040 == 0, 'entry IEEE32 execution state mismatch'",
     "            if not self.is_preflight:\n                assert row['mxcsr'] & 0xe040 == 0, 'entry IEEE32 execution state mismatch'"),
    ("        breakpoints.append(Entered('sirius_native_' + slug + '_try_packet'))",
     "        breakpoints.append(Entered('sirius_native_' + slug + '_try_packet'))\n        breakpoints.append(Entered('sirius_native_' + slug + '_admit_packet', True))"),
]
for before, after in changes:
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
