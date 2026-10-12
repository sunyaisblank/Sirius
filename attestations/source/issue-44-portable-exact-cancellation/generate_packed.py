"""Independent lattice rounding for the original ordered pack contracts."""
from pathlib import Path
from collections import Counter
import hashlib
import importlib.util
import itertools
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
PRIVATE_BOUND = 0x7d800000
MASK = ref.MASK
counts = Counter()
records = []

def append(op, a, b, c, error, family):
    if op == 2:
        assert all((v & MASK) <= PRIVATE_BOUND for v in (a, b, c))
    records.append((op, a, b, c, error, 0, 0, 0))
    counts[family] += 1

edges = [v | sign for v in (0, 1, 2, 0x7fffff, 0x800000, 0x0bffffff,
          0x0c000000, 0x33800000, 0x3f800000, 0x3f800001, BOUND - 1, BOUND,
          BOUND + 1, ref.INF, ref.INF + 1, ref.NAN) for sign in (0, ref.SIGN)]
for a, b in itertools.product(edges, repeat=2):
    append(0, a, b, 0, 0, 'pair_boundary_cross')
for sign in (0, ref.SIGN):
    for exponent in range(1, 248):
        for mantissa in (0, 1, 0x400000, 0x7fffff):
            a = ((exponent << 23) | mantissa) | sign
            for delta in (-1, 0, 1):
                b = ((a & MASK) + delta) | (sign ^ ref.SIGN)
                append(0, a, b, 0, 0, 'pair_opposite_and_adjacent')
                append(1, a, b, 1, 0, 'triple_opposite_and_adjacent')
for limbs in itertools.product((0, ref.SIGN, 1, 0x800000, 0x3f800000, 0xbf800000), repeat=3):
    append(1, *limbs, 0, 'triple_zero_sign_and_cancellation_cross')
for op in (0, 1, 2):
    for error in (0, 1, 0x800000, BOUND, BOUND + 1, ref.SIGN, 0x80000001, ref.INF, ref.NAN):
        for a in (0, 0x3f800000, BOUND):
            append(op, a, a ^ ref.SIGN, 1, error, 'radius_refusal_and_boundary')
for sign in (0, ref.SIGN):
    for a in (BOUND, BOUND + 1, PRIVATE_BOUND - 1, PRIVATE_BOUND):
        for tail in (0, 1, 0x3f800000, BOUND):
            append(2, a | sign, a | (sign ^ ref.SIGN), tail, 0, 'private_large_term_cancellation')
rng = random.Random(0x5349524955530075)
for _ in range(2048):
    limbs = [rng.randrange(BOUND + 1) | rng.choice((0, ref.SIGN)) for _ in range(3)]
    append(1, *limbs, rng.choice((0, 1, 0x800000)), 'bounded_raw_three_limb_fixed_seed')

def sum32(a, b, mode):
    # The public word contract fixes magnitude ordering and three separate RTE
    # operations. Each operation uses exact Fraction lattice rounding, without
    # the production significand/shift-jam algorithm or native floating point.
    if (a & MASK) < (b & MASK):
        a, b = b, a
    am, bm = a & MASK, b & MASK
    if mode == 'normal_transport' and bm == 0 and am < ref.INF:
        return (a if am else a & b), b
    high = ref.add(a, b)
    displaced = ref.add(high, a ^ ref.SIGN)
    low = ref.add(b, displaced ^ ref.SIGN)
    if mode != 'integer' and 0x0c000000 <= am <= 0x7e800000 and 0x0c000000 <= bm <= 0x7e800000:
        high = 0 if high & MASK == 0 else high
        low = 0 if low & MASK == 0 else low
    return high, low

def expected(row, mode):
    op, a, b, c, error, *_ = row
    invalid = [0] * 8
    if (error & MASK) > BOUND or error & ref.SIGN:
        return invalid
    if op == 0:
        if (a & MASK) > BOUND or (b & MASK) > BOUND:
            return invalid
        high, low = sum32(a, b, mode)
        tail = 0
    else:
        if op == 1 and any((v & MASK) > BOUND for v in (a, b, c)):
            return invalid
        # Original five magnitude-ordered error-free transforms, with distinct
        # intermediate rounding/sign words. Their final exact real sum is also
        # checked independently below, rather than inferred from word agreement.
        first, first_residual = sum32(b, c, mode)
        second, second_residual = sum32(a, first, mode)
        third, third_residual = sum32(second_residual, first_residual, mode)
        high, carry = sum32(second, third, mode)
        low, tail = sum32(carry, third_residual, mode)
    if (high & MASK) > BOUND:
        return invalid
    total = ref.exact(a) + ref.exact(b) + (ref.exact(c) if op else 0)
    assert ref.exact(high) + ref.exact(low) + ref.exact(tail) == total
    return [high, low, tail, error, 1, 0, 0, 0]

assert len(records) <= 131071
inputs = [len(records)] + [v for row in records for v in row]
assert 4 * (len(inputs) + 8 * len(records)) <= 8 * 1024 * 1024
(work / 'packed-input.bin').write_bytes(struct.pack('<%dI' % len(inputs), *inputs))
identity = lambda p: {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
report = {'records': len(records), 'partitions': dict(counts), 'input': identity(work / 'packed-input.bin'),
          'complete_output_words_each': 8 * len(records), 'requested_device_bytes': 4 * (len(inputs) + 8 * len(records)),
          'modes': {}, 'expectation_method': 'Unchanged independent Fraction lattice rounding of the pre-existing ordered pack contract; every valid result additionally preserves the exact real limb sum and unchanged radius. Complete status and padding words included.',
          'scope': 'Finite RPPack/RTPack/private RTFinish witnesses. No full pipeline, controller, benefit or qualification claim.'}
for mode in ('integer', 'normal', 'normal_transport'):
    values = [expected(row, mode) for row in records]
    path = work / ('packed-expected-' + mode + '.bin')
    flat = [v for value in values for v in value]
    path.write_bytes(struct.pack('<%dI' % len(flat), *flat))
    report['modes'][mode] = {'expected': identity(path), 'valid_records': sum(v[4] for v in values),
                             'complete_zero_refusals': sum(v[4] == 0 for v in values)}
(work / 'packed-fixture.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report))
