"""Input partitions for word preservation plus independent quotient enclosures."""
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
counts = Counter()
cases = []


def bits(e):
    assert -149 <= e <= 120
    return (e + 127) << 23 if e >= -126 else 1 << (e + 149)


def append(a, b, family, status=-1, contract='interval'):
    assert len(a) == len(b) == 4 and all(0 <= v <= 0xffffffff for v in (*a, *b))
    cases.append({'a': list(a), 'b': list(b), 'family': family,
                  'required_status': status, 'contract': contract})
    counts[family] += 1


unit = (0x3f800000, 0, 0, 0)
recovery_a, recovery_b = unit, (0x3fc00000, 0, 0, 0)
# Scalar self-division, unit-denominator identity and zero numerators are
# supported algebraic references, independent of the implementation's words.
for e in (-149, -148, -127, -126, -103, -60, -30, -1, 0, 1, 30, 60):
    for sign in (0, ref.SIGN):
        x = bits(e) | sign
        append((x, 0, 0, 0), (x, 0, 0, 0), 'scalar_self_division', 1, 'self')
        append((x, 0, 0, 0), unit, 'unit_denominator_scalar_identity', 1, 'unit')
        append((0, ref.SIGN, 0, 0), (x, 0, 0, 0), 'signed_zero_numerator', 1, 'zero')
for signs in itertools.product((0, ref.SIGN), repeat=3):
    append((*signs, 0), recovery_b, 'all_zero_limb_signs', 1, 'zero')

for ea in (-60, -30, -10, -1, 0, 1, 10, 30, 60):
    for eb in (-30, -10, 0, 10, 30):
        for sign in (0, ref.SIGN):
            for offset in (-1, 0, 1):
                a = (bits(ea) + offset) | sign
                b = bits(eb) + 1
                append((a, bits(max(-149, ea - 25)), bits(max(-149, ea - 50)), 0),
                       (b, bits(max(-149, eb - 25)), bits(max(-149, eb - 50)), 0),
                       'nonzero_three_limb_corrections_and_adjacent_binades')
                append((a, bits(max(-149, ea - 25)) | ref.SIGN,
                        bits(max(-149, ea - 50)), bits(max(-149, ea - 55))),
                       (b | ref.SIGN, bits(max(-149, eb - 25)),
                        bits(max(-149, eb - 50)) | ref.SIGN, bits(max(-149, eb - 55))),
                       'signed_three_limb_uncertainty')
    for sign in (0, ref.SIGN):
        append((bits(ea) | sign, bits(max(-149, ea - 25)),
                bits(max(-149, ea - 50)) | ref.SIGN, bits(max(-149, ea - 55))),
               unit, 'unit_denominator_full_pack_identity', 1, 'unit')

for e in (-149, -148, -127, -126, -103, -60, -1, 0, 1, 60, 120):
    for f in (-152, -151, -150, -149, -127, -126, -125, -1, 0, 1, 119, 120, 121):
        ae = e + f
        if -149 <= ae <= 120:
            for sa in (0, ref.SIGN):
                for sb in (0, ref.SIGN):
                    append((bits(ae) | sa, 0, 0, 0), (bits(e) | sb, 0, 0, 0),
                           'complete_division_seed_exponent_and_subnormal_boundaries')

bad = (BOUND + 1, ref.INF, ref.NAN, ref.SIGN, 0x80000001)
for value in bad:
    for component in range(4):
        a, b = list(unit), list(unit)
        # Sign words in limb positions are valid zeros; only use them as
        # negative radii. Nonfinite/out-of-bound limb words always refuse.
        if value in (ref.SIGN, 0x80000001) and component != 3:
            continue
        a[component] = value
        append(a, recovery_b, 'invalid_numerator_pack', 0, 'refusal')
        b[component] = value
        append(recovery_a, b, 'invalid_denominator_pack', 0, 'refusal')
for center in (0, 1, 0x00800000, 0x3f800000):
    for sign in (0, ref.SIGN):
        append(recovery_a, (center | sign, 0, 0, center),
               'denominator_interval_contains_zero', 0, 'refusal')
for b in ((0x3f800000, 0xbf800000, 0, 0),
          (0x3f800000, 0xbf800000, 1, 1),
          (0xbf800000, 0x3f800000, ref.SIGN, 0)):
    append(recovery_a, b, 'denominator_cancellation_interval', 0, 'refusal')

rng = random.Random(0x5349524955530077)
for _ in range(4096):
    ea, eb = rng.randrange(-90, 41), rng.randrange(-40, 41)
    a = ((bits(ea) | rng.randrange(1 << 23)) | rng.choice((0, ref.SIGN)),
         bits(max(-149, ea - 25)) | rng.choice((0, ref.SIGN)),
         bits(max(-149, ea - 50)) | rng.choice((0, ref.SIGN)),
         rng.choice((0, 1, bits(max(-149, ea - 60)))))
    b = ((bits(eb) | rng.randrange(1 << 23)) | rng.choice((0, ref.SIGN)),
         bits(max(-149, eb - 25)) | rng.choice((0, ref.SIGN)),
         bits(max(-149, eb - 50)) | rng.choice((0, ref.SIGN)),
         rng.choice((0, 1, bits(max(-149, eb - 60)))))
    append(a, b, 'fixed_seed_canonical_scale_raw_triples')

# Adjacent refusal/recovery and cross-workgroup repeated values expose any
# accidental preparation state escaping one operation. Every marker row is
# the same exact rational2/3 input and must reproduce all complete words.
original = list(cases)
cases = []
recovery_indices = []
for i, case in enumerate(original):
    if i % 31 == 0:
        recovery_indices.append(len(cases))
        append(recovery_a, recovery_b, 'recovery_marker', 1, 'recovery')
    cases.append(case)
    if case['required_status'] == 0:
        recovery_indices.append(len(cases))
        append(recovery_a, recovery_b, 'recovery_marker', 1, 'recovery')
recovery_indices.append(len(cases))
append(recovery_a, recovery_b, 'recovery_marker', 1, 'recovery')
inputs = [len(cases)]
for case in cases:
    inputs.extend((0, *case['a'], *case['b'], *([0] * 7)))
words = 8 * len(cases)
requested = 4 * (len(inputs) + words)
assert len(cases) <= 87381 and requested <= 8 * 1024 * 1024
for name, values in (('triple-input', inputs), ('triple-capture-shape', [0] * words)):
    p = work / (name + '.bin')
    assert not p.exists()
    p.write_bytes(struct.pack('<%dI' % len(values), *values))
identity = lambda p: {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
report = {'records': len(cases), 'complete_words_each': words, 'partitions': dict(counts),
          'requested_device_bytes': requested, 'recovery_indices': recovery_indices,
          'input': identity(work / 'triple-input.bin'), 'cases': cases,
          'expectation_method': 'Original/candidate complete-word preservation plus independent exact input quotient-interval enclosure, raw bounds/refusals, unit pack identity and zero/self/recovery obligations. Captured original words are a comparison target, not an independent numerical oracle.',
          'scope': 'Finite complete retained division/refusal/recovery witnesses; integrated scientific/stream contracts, benefit and qualification remain separate.'}
(work / 'triple-fixture.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({k: v for k, v in report.items() if k not in ('cases', 'recovery_indices')}))
