"""Independent seed expectations retain the two distinct scaling rounds."""
from pathlib import Path
from collections import Counter
from fractions import Fraction
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
BOUND = 0x7b800000
counts, outcomes = Counter(), Counter()
rows, expected = [], []


def power(e):
    return Fraction(1 << e) if e >= 0 else Fraction(1, 1 << -e)


def seed_expected(a, b):
    if a & ref.MASK == 0:
        outcomes['zero_numerator_positive_zero_precedence'] += 1
        return 0
    aa, bb = abs(ref.exact(a)), abs(ref.exact(b))
    ea, eb = ref.exponent_of(aa), ref.exponent_of(bb)
    exponent = ea - eb
    if exponent < -151 or exponent > 120:
        outcomes['original_exponent_refusal_positive_zero'] += 1
        return 0
    normalized_ratio = (aa / power(ea)) / (bb / power(eb))
    seed = ref.round_exact(normalized_ratio)
    if exponent < -126:
        outcomes['both_distinct_scaling_rounds'] += 1
        seed = ref.round_exact(ref.exact(seed) * power(-126))
        exponent += 126
    seed = ref.round_exact(ref.exact(seed) * power(exponent))
    result = seed ^ ((a ^ b) & ref.SIGN)
    if result & ref.MASK == 0:
        outcomes['accepted_negative_zero' if result else 'accepted_positive_zero'] += 1
    return result


def append(a, b, family):
    assert a & ref.MASK <= BOUND and 0 < b & ref.MASK <= BOUND
    rows.append((0, a, b, 0, 0, 0, 0, 0))
    expected.extend((seed_expected(a, b), 0, 0, 0, 1, 0, 0, 0))
    counts[family] += 1


subnormals = {v for bit in range(23) for v in ((1 << bit) - 1, 1 << bit, (1 << bit) + 1)
              if 0 < v < 0x00800000}
denominators = sorted(subnormals | {((e << 23) | m)
                       for e in range(1, 248) for m in (0, 1, 0x7fffff)
                       if ((e << 23) | m) <= BOUND})
differences = (-152, -151, -150, -149, -128, -127, -126, -125, -1, 0, 1, 119, 120, 121)
for b in denominators:
    eb = ref.exponent_of(ref.exact(b))
    for difference in differences:
        ea = eb + difference
        if not -149 <= ea <= 120:
            continue
        # Input selection is raw boundary/adjacency enumeration; all expected
        # results above are derived from independently decoded exact values.
        boundary = (ea + 127) << 23 if ea >= -126 else 1 << (ea + 149)
        for delta in (-1, 0, 1):
            a = boundary + delta
            if not 0 < a <= BOUND:
                continue
            for sa in (0, ref.SIGN):
                for sb in (0, ref.SIGN):
                    append(a | sa, b | sb, 'all_denominator_exponents_subnormal_bits_and_adjacent_seed_cutoffs')
    for sa in (0, ref.SIGN):
        for sb in (0, ref.SIGN):
            append(sa, b | sb, 'signed_zero_numerator_cross')

rng = random.Random(0x5349524955530076)
for _ in range(16384):
    a = rng.randrange(BOUND + 1) | rng.choice((0, ref.SIGN))
    b = rng.randrange(1, BOUND + 1) | rng.choice((0, ref.SIGN))
    append(a, b, 'fixed_seed_bounded_finite_pairs')

inputs = [len(rows)] + [v for row in rows for v in row]
requested = 4 * (len(inputs) + len(expected))
assert len(rows) <= 131071 and requested <= 8 * 1024 * 1024
identity = lambda p: {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
for name, words in (('seed-input', inputs), ('seed-expected', expected)):
    p = work / (name + '.bin')
    assert not p.exists()
    p.write_bytes(struct.pack('<%dI' % len(words), *words))
record = {'records': len(rows), 'complete_words_each': len(expected), 'partitions': dict(counts),
          'outcomes': dict(outcomes), 'requested_device_bytes': requested,
          'input': identity(work / 'seed-input.bin'), 'expected': identity(work / 'seed-expected.bin'),
          'expectation_method': 'Exact Fraction mantissa quotient independently rounded, followed by each original scaling round separately and original zero/exponent/sign policy. No production word arithmetic or native float generates expectations.',
          'scope': 'Finite prepared seed contract only; complete triple division/stream/refusal/recovery, benefit and qualification remain separate.'}
(work / 'seed-fixture.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record))
