"""Keep the original oracle prefix and independently round normalized quotients."""
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
minimum, maximum = 1 << 23, (1 << 24) - 1


def add(a, d, family):
    if not minimum <= a <= maximum or not minimum <= d <= maximum:
        return
    aw, dw = 0x3f800000 | (a - minimum), 0x3f800000 | (d - minimum)
    # Exact rational quotient rounded to the binary32 lattice, independently
    # of the production reciprocal, chunks, guard/sticky bits and Round.
    high = ref.round_exact(ref.exact(aw) / ref.exact(dw))
    rows.append((14, aw, dw, high, 0, 1))
    counts[family] += 1


edges = [minimum, minimum + 1, minimum + 2, 0xbfffff, 0xc00000,
         0xc00001, maximum - 2, maximum - 1, maximum]
for a in edges:
    for d in edges:
        add(a, d, 'mantissa_boundary_cross')

denominators = set(edges)
for reciprocal in range(256, 512):
    edge = 0xffffffff // reciprocal
    for delta in (-1, 0, 1):
        d = edge + delta
        if minimum <= d <= maximum:
            denominators.add(d)
for d in sorted(denominators):
    for a in [*edges, d - 1, d, d + 1]:
        add(a, d, 'reciprocal_bucket_boundaries_and_normalization_adjacent')
    # Test numerators adjacent to exact output-lattice midpoints. Positive
    # normalized input quotients cannot be exact halfway cases: after gcd,
    # a dyadic denominator has at most23 factors of two, while a midpoint
    # in these output binades needs24 or25. Near-midpoints remain relevant.
    for shift in (24, 25):
        for q in (minimum, minimum + 1, 0xbfffff, 0xc00000, maximum - 1, maximum):
            a = d * (2 * q + 1) // (1 << shift)
            for delta in (-1, 0, 1):
                add(a + delta, d, 'adjacent_to_exact_output_midpoints')

for bit in range(23):
    for delta in (-1, 0, 1):
        d = minimum + (1 << bit) + delta
        for a in [minimum, maximum, d - 1, d, d + 1]:
            add(a, d, 'fraction_bit_boundaries')

rng = random.Random(0x5349524955530075)
for _ in range(32768):
    add(rng.randrange(minimum, maximum + 1), rng.randrange(minimum, maximum + 1),
        'fixed_seed_normalized_pairs')

inputs = [len(rows)]
expected = []
for op, a, b, high, low, valid in rows:
    inputs.extend((op, a, b, 0))
    expected.extend((high, low, valid, 0))
requested = 4 * (len(inputs) + len(expected))
assert len(rows) <= 262143 and requested <= 8 * 1024 * 1024
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
          'expectation_method': 'Unchanged independent exact Fraction-lattice oracle; complete original prefix and order preserved. Operation14 has positive normalized inputs only. No native float or production word algorithm generates expectations.',
          'scope': 'Finite normalized primitive partitions; prepared seed scaling/complete retained divisions, integrated controls, benefit and qualification remain separate.'}
(work / 'fixture.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record))
