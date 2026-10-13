"""Project the ordered independent residual corpus without model injection."""
from pathlib import Path
import hashlib
import importlib.util
import json
import resource
import struct
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024**2, 64 * 1024**2))
root = Path.cwd().resolve()
work = Path(__file__).resolve().parent
reference_path = root / 'tests/support/portable_binary32/reference.py'
fixture_path = root / 'bin/linux-gcc/tests/backend/portable_binary32_reference.bin'
assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == '6b6474bf66c86713980cba81f0cee1b290f54c8f42700ca9353aa108cf0c185e'
raw = fixture_path.read_bytes()
assert hashlib.sha256(raw).hexdigest() == '3fb1433f9cd14fbc9fd2b99a38a4893a4a571150df10f2fbe58172318fc7691e'
spec = importlib.util.spec_from_file_location('independent_residual_reference', reference_path)
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
ref.selfcheck()
magic, version, count = struct.unpack_from('<III', raw)
assert (magic, version, count) == (ref.MAGIC, ref.VERSION, ref.RECORD_COUNT)
assert len(raw) == 12 + 24 * count
records = [row for row in struct.iter_unpack('<6I', raw[12:]) if row[0] == 13]
assert len(records) == 53524
assert [(row[1], row[2]) for row in records] == [(a, b) for a, b, _ in ref.residual_cases()]
witnesses = []
for sa in (0, 0x80000000):
    for sb in (0, 0x80000000):
        a, b = 0x3f7ff400 | sa, 0x3f800600 | sb
        high, low, valid = ref.product_expected(a, b)
        witnesses.append((13, a, b, high, low, valid))
records += witnesses
inputs = [len(records)]
expected = []
for op, a, b, high, low, valid in records:
    assert op == 13
    inputs += [op, a, b, 0]
    expected += [high, low, valid, 0]
projected = [('residual-input.bin', struct.pack('<%dI' % len(inputs), *inputs),
              '4573b5bf349f1510652f75361528194184c74a61ad8dd0280273b5de036f2409'),
             ('residual-expected.bin', struct.pack('<%dI' % len(expected), *expected),
              '2110fe70214d4325c01bd36a4f27ff95acbf4f34e1816f05a31508d22dcce497')]
identities = {}
for name, data, digest in projected:
    assert hashlib.sha256(data).hexdigest() == digest
    path = work / name
    assert not path.exists()
    path.write_bytes(data)
    identities[name] = {'bytes': len(data), 'sha256': digest}
report = {'records': len(records), 'original_ordered_records': 53524,
    'additional_signed_binade_witnesses': witnesses, 'complete_output_words': len(expected),
    'requested_device_bytes': 4 * (len(inputs) + len(expected)),
    'groups': (len(records) + 63) // 64, 'files': identities,
    'input_fields': ['opcode13', 'a', 'b', 'zero-padding'],
    'output_fields': ['high', 'low', 'valid', 'zero-padding'],
    'hardware_estimate_input_present': False,
    'scope': 'Ordered independent residual expectation projection plus four signed witnesses. Shader execution remains required; no injected hardware estimate.'}
(work / 'residual-fixture.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report))
