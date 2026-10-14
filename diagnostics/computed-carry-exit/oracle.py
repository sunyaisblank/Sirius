"""Deterministic Fraction expectations for the original five-word accumulator.

No candidate implementation is evaluated. Finite TwoSum residuals are derived
from exact real sums and independently cross-checked against three rounded adds.
"""
from collections import Counter
from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import random
import struct
import sys

sys.dont_write_bytecode = True
ROOT = Path.cwd().resolve()
WORK = Path(__file__).resolve().parent
REFERENCE = ROOT / 'tests/support/portable_binary32/reference.py'
spec = importlib.util.spec_from_file_location('issue101_fraction_reference', REFERENCE)
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
SIGN, MASK, PUBLIC, PRIVATE = 0x80000000, 0x7fffffff, 0x7b800000, 0x7d800000
finite_splits = 0


def two_sum(a, b):
    global finite_splits
    if (a & MASK) < (b & MASK):
        a, b = b, a
    high = ref.add(a, b)
    displacement = ref.add(high, a ^ SIGN)
    chain_low = ref.add(b, displacement ^ SIGN)
    if ref.exact(a) is not None and ref.exact(b) is not None and ref.exact(high) is not None:
        low = b if (b & MASK) == 0 else ref.round_exact(ref.exact(a) + ref.exact(b) - ref.exact(high))
        assert low == chain_low, (a, b, high, low, chain_low)
        finite_splits += 1
    else:
        low = chain_low  # Exceptional intermediate; final private guards own refusal.
    return high, low


def upward(a, b):
    if (a & MASK) > PRIVATE or (b & MASK) > PRIVATE:
        return 0x7f000000
    if (a & MASK) == 0:
        return b
    if (b & MASK) == 0:
        return a
    return (ref.add(a, b) + 1) & 0xffffffff


def expected(a, x):
    if a[4] == 0 or (x & MASK) > PUBLIC:
        return (0, 0, 0, 0, 0), []
    carry = x
    limbs, carries = [], []
    for limb in a[:3]:
        upper, carry = two_sum(limb, carry)
        limbs.append(upper)
        carries.append(carry)
    radius = upward(a[3], carry & MASK)
    if any((v & MASK) > PRIVATE for v in (*limbs, radius)):
        return (0, 0, 0, 0, 0), carries
    return (*limbs, radius, 1), carries


def power(e):
    return ref.round_exact(Fraction(2**e) if e >= 0 else Fraction(1, 2**-e))


def records():
    rows = []
    def add(label, a, x):
        assert len(a) == 5
        rows.append((label, tuple(a), x))
    first = ((0, 0, 0, 0, 1), 0x3f800000)
    second = ((0x3f800000, power(-24) ^ SIGN, 0, 0, 1), power(-24))
    full = ((0x3f800000, 0x3f800000, 0x3f800000, 1, 1), power(-24))
    zero = ((SIGN, 0, SIGN, SIGN, 1), SIGN)
    # Adjacent eight-record cohorts alternate the new cuts and unchanged paths.
    for _ in range(16):
        for a, x in (first, full, second, zero, full, first, zero, second):
            add('mixed_adjacent_records', a, x)
    remaining = (0, SIGN, 1, SIGN | 1, 0x007fffff, 0x807fffff,
                 0x3f800000, 0xbf800000, PRIVATE, PRIVATE | SIGN,
                 PRIVATE + 1, (PRIVATE + 1) | SIGN,
                 0x7f800000, 0xff800000, 0x7f800001, 0xffc00001)
    initial = (1, SIGN | 1, 0x00800000, 0x80800000,
               0x3f800000, 0xbf800000, PUBLIC, PUBLIC | SIGN)
    for high, x, low, tail in itertools.product((0, SIGN), initial, remaining, remaining):
        add('first_computed_zero', (high, low, tail, 0, 1), x)
    for e, gap, hs, xs in itertools.product((-100, -20, 0, 60, 100), (23, 24, 25, 26), (0, SIGN), (0, SIGN)):
        high, x = power(e) | hs, power(e - gap) | xs
        _, residual = two_sum(high, x)
        if residual & MASK:
            for tail in remaining:
                add('second_computed_zero', (high, residual ^ SIGN, tail, 0, 1), x)
            add('nonzero_carry_through_three', (high, high, high, 1, 1), x)
    for high, low, tail, radius, x in itertools.product(
            (0, SIGN, 0x3f800000, 0xbf800000), (0, SIGN), (0, SIGN), (0, SIGN, 1), (0, SIGN)):
        add('initial_zero_unchanged', (high, low, tail, radius, 1), x)
    centers = (0, SIGN, PRIVATE - 1, PRIVATE, PRIVATE + 1,
               (PRIVATE - 1) | SIGN, PRIVATE | SIGN, (PRIVATE + 1) | SIGN,
               0x7f7fffff, 0xff7fffff, 0x7f800000, 0xff800000, 0x7fc00000, 0xff800001)
    radii = (0, SIGN, 1, 0x007fffff, 0x00800000, PRIVATE - 1, PRIVATE, PRIVATE + 1,
             0xbf800000, PRIVATE | SIGN, 0x7f800000, 0x7fc00000)
    for position, limb, radius in itertools.product(range(3), centers, radii):
        a = [0, 0, 0, radius, 1]
        a[position] = limb
        for x in (0, SIGN, 1, 0x3f800000, PUBLIC):
            add('private_and_radius_boundaries', a, x)
    for valid, x, high in itertools.product((0, 1, 0xffffffff),
            (0, SIGN, PUBLIC - 1, PUBLIC, PUBLIC + 1, (PUBLIC + 1) | SIGN,
             0x7f800000, 0xff800000, 0x7f800001, 0xffc00001),
            (0, SIGN, 0x3f800000, 0x7f800000, 0x7fc00000)):
        add('initial_validity_scalar_refusal', (high, 0, SIGN, 0, valid), x)
    # Public-input ties and cancellation with overlapping private limbs.
    edges = (1, 2, 3, 0x003fffff, 0x007fffff, 0x00800000, 0x00800001,
             0x33800000, 0x3f7fffff, 0x3f800000, 0x3f800001, PUBLIC - 1, PUBLIC)
    for h, x, sign in itertools.product(edges, edges, (0, SIGN)):
        add('subnormal_tie_cancellation', (h, h ^ SIGN, 1, 1, 1), x | sign)
    rng = random.Random(10120261014)
    for _ in range(2048):
        limbs = [(rng.randrange(PRIVATE + 1) | rng.choice((0, SIGN))) for _ in range(3)]
        add('fixed_seed_finite_overlap', (*limbs, rng.randrange(PRIVATE + 1), 1),
            rng.randrange(PUBLIC + 1) | rng.choice((0, SIGN)))
    return rows


def main():
    ref.selfcheck()
    assert expected((0, 0, 0, 0, 1), 0x3f800000)[0] == (0x3f800000, 0, 0, 0, 1)
    assert expected((SIGN, SIGN, SIGN, SIGN, 1), 0x3f800000)[0] == (0x3f800000, SIGN, SIGN, 0, 1)
    assert expected((0, 0, 0, 0, 0), 0x3f800000)[0] == (0, 0, 0, 0, 0)
    rows = records()
    inputs, outputs = [len(rows)], []
    partitions, cuts, signs = Counter(), Counter(), Counter()
    for label, a, x in rows:
        words, carries = expected(a, x)
        inputs.extend((*a, x))
        outputs.extend(words)
        partitions[label] += 1
        if (x & MASK) != 0 and carries:
            if (carries[0] & MASK) == 0:
                cuts['after_first_sum'] += 1
                signs['first_negative_zero' if carries[0] & SIGN else 'first_positive_zero'] += 1
            elif (carries[1] & MASK) == 0:
                cuts['after_second_sum'] += 1
            else:
                cuts['no_early_zero'] += 1
        if words[4] == 0:
            assert words == (0, 0, 0, 0, 0)
    assert cuts['after_first_sum'] and cuts['after_second_sum'] and cuts['no_early_zero']
    assert signs['first_negative_zero'] and signs['first_positive_zero']
    raw_input = struct.pack('<%dI' % len(inputs), *inputs)
    raw_expected = struct.pack('<%dI' % len(outputs), *outputs)
    assert len(raw_input) + len(raw_expected) < 8 * 1024**2
    (WORK / 'input.bin').write_bytes(raw_input)
    (WORK / 'expected.bin').write_bytes(raw_expected)
    meta = {'records': len(rows), 'input_words_per_record': 6, 'output_words_per_record': 5,
        'input_header_words': 1, 'complete_output_words': len(outputs),
        'requested_device_bytes': len(raw_input) + len(raw_expected),
        'partitions': dict(partitions), 'original_computed_carry_paths': dict(cuts),
        'first_computed_zero_signs': dict(signs), 'finite_exact_splits_crosschecked': finite_splits,
        'reference_sha256': hashlib.sha256(REFERENCE.read_bytes()).hexdigest(),
        'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'input_sha256': hashlib.sha256(raw_input).hexdigest(),
        'expected_sha256': hashlib.sha256(raw_expected).hexdigest(),
        'expected_modes': ['integer', 'normal-transport'],
        'scope': 'Finite independent Fraction oracle for the original accumulator, shared by both frozen producers/modes. No candidate execution, compiler, provider, integrated-control or benefit verdict.'}
    (WORK / 'fixture.json').write_text(json.dumps(meta, indent=2) + '\n')
    print(json.dumps(meta, indent=2))


if __name__ == '__main__':
    main()
