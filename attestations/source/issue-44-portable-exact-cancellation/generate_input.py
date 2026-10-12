"""Preserve the full original oracle and add exact cancellation partitions."""
from pathlib import Path
from collections import Counter
import hashlib
import importlib.util
import json
import random
import struct
import sys

sys.dont_write_bytecode = True
root = Path.cwd()
work = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('independent_reference', root / 'tests/support/portable_binary32/reference.py')
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
ref.selfcheck()
fixture = root / 'bin/linux-gcc/tests/backend/portable_binary32_reference.bin'
raw = fixture.read_bytes()
assert struct.unpack_from('<3I', raw) == (0x32334250, 1, 174268)
assert len(raw) == 12 + 24 * 174268
rows = list(struct.iter_unpack('<6I', raw[12:]))
counts = Counter()

def add(op, a, b, family):
    rows.append((op, a, b, ref.scalar_expected(op, a, b), 0, 1))
    counts[family] += 1

normal = [(e << 23) | m for e in range(1, 255)
          for m in (0, 1, 2, 0x3fffff, 0x400000, 0x400001, 0x7ffffe, 0x7fffff)]
subnormal = sorted({n for bit in range(23) for n in ((1 << bit) - 1, 1 << bit, (1 << bit) + 1)
                    if 0 < n < 0x800000} | {0x7ffffe, 0x7fffff})
for family, magnitudes in (('all_normal_exponents_and_adjacent', normal), ('subnormal_leading_bits_and_adjacent', subnormal)):
    for magnitude in magnitudes:
        for sign in (0, ref.SIGN):
            a = magnitude | sign
            for delta in (-1, 0, 1):
                neighbor = magnitude + delta
                assert 0 <= neighbor <= ref.INF
                add(0, a, neighbor | (sign ^ ref.SIGN), family)
                add(1, a, neighbor | sign, family)

special = [v | sign for v in (0, ref.INF, ref.INF + 1, ref.NAN, ref.NAN + 1, 0x7fffffff)
           for sign in (0, ref.SIGN)]
for a in special:
    for b in special:
        add(0, a, b, 'zero_infinity_nan_precedence')
        add(1, a, b, 'zero_infinity_nan_precedence')

rng = random.Random(0x5349524955530074)
for _ in range(8192):
    magnitude = rng.randrange(1, ref.INF)
    for sign in (0, ref.SIGN):
        a = magnitude | sign
        add(0, a, a ^ ref.SIGN, 'fixed_seed_finite_opposite_equal')
        add(1, a, a, 'fixed_seed_finite_opposite_equal')

assert len(rows) <= 262143
inputs = [len(rows)]
expected = []
for op, a, b, high, low, valid in rows:
    inputs.extend((op, a, b, 0))
    expected.extend((high, low, valid, 0))
requested = 4 * (len(inputs) + len(expected))
assert requested <= 8 * 1024 * 1024
identity = lambda p: {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
for name, words in (('input', inputs), ('expected', expected)):
    target = work / (name + '.bin')
    assert not target.exists()
    target.write_bytes(struct.pack('<%dI' % len(words), *words))
record = {'original_records': 174268, 'additional_records': len(rows) - 174268,
          'partitions': dict(counts), 'records': len(rows), 'complete_expected_words': len(expected),
          'requested_device_bytes': requested, 'device_allocation_limit': 8 * 1024 * 1024,
          'original_fixture': identity(fixture), 'input': identity(work / 'input.bin'),
          'expected': identity(work / 'expected.bin'),
          'expectation_method': 'Unchanged independent exact Fraction lattice oracle; full original record sequence retained as the exact prefix. No native float arithmetic or production word algorithm generates expectations.',
          'scope': 'Finite primitive partitions only. Connected complete packed sums, products, integrated streams, benefit and qualification remain separate.'}
(work / 'fixture.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record))
