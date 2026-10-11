"""Screen a conditional native-default TwoProduct substitution, without GPU work.

Fraction arithmetic models each original ordered binary32 operation. The existing
independent corpus supplies exact product/residual expectations. Matching this
finite set is not exhaustive equivalence or actual device/compiler qualification.
"""
from pathlib import Path
from collections import Counter
import hashlib
import importlib.util
import json
import struct
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
MASK, SIGN = 0x7fffffff, 0x80000000
LOWER, UPPER = 0x2b800000, 0x5d800000
REFERENCE = ROOT / 'tests/support/portable_binary32/reference.py'
CORPUS = ROOT / 'bin/linux-gcc/tests/backend/portable_binary32_reference.bin'
SOURCE = ROOT / 'src/sirius/kernels/retained_pair.slang'
spec = importlib.util.spec_from_file_location('independent_binary32_reference', REFERENCE)
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)

def seal(path):
    raw = path.read_bytes()
    return {'path': str(path.relative_to(ROOT)), 'bytes': len(raw),
            'sha256': hashlib.sha256(raw).hexdigest()}

def ordered(op, a, b):
    return ref.scalar_expected(op, a, b)

def split(a):
    c = ordered(2, 0x45800800, a)
    high = ordered(1, c, ordered(1, c, a))
    return high, ordered(1, a, high)

def original(a, b):
    ah, al = split(a)
    bh, bl = split(b)
    high = ordered(2, a, b)
    e1 = ordered(1, high, ordered(2, ah, bh))
    e2 = ordered(1, e1, ordered(2, al, bh))
    e3 = ordered(1, e2, ordered(2, ah, bl))
    return high, ordered(1, ordered(2, al, bl), e3)

def fused_residual(a, b, high):
    # A genuinely fused instruction rounds only the exact final sum. All
    # admitted operands/results are finite; any zero sign is rejected below.
    return ref.round_exact(ref.exact(a) * ref.exact(b) - ref.exact(high))

def candidate(a, b):
    if LOWER <= (a & MASK) <= UPPER and LOWER <= (b & MASK) <= UPPER:
        high = ordered(2, a, b)
        low = fused_residual(a, b, high)
        if low & MASK:
            return (high, low), 'guarded-nonzero'
        return original(a, b), 'zero-residual-fallback'
    return original(a, b), 'outside-guard-fallback'

def check(a, b, label, expected=None):
    if expected is not None:
        assert expected == ref.product_expected(a, b), 'corpus/oracle discrepancy'
    if (a & MASK) > UPPER or (b & MASK) > UPPER:
        assert expected is None or expected == (0, 0, 0)
        counts['caller-refused'] += 1
        return
    old = original(a, b)
    new, route = candidate(a, b)
    counts[route] += 1
    families[label] += 1
    assert new == old, (label, hex(a), hex(b), old, new, route)
    exact_product = ref.exact(a) * ref.exact(b)
    assert old[0] == ref.round_exact(exact_product, (a ^ b) & SIGN)
    error = abs(exact_product - ref.exact(old[0]) - ref.exact(old[1]))
    assert error <= 5 * ref.exact(1), (label, hex(a), hex(b), error)
    if route == 'guarded-nonzero':
        assert new[1] & MASK >= 0x00800000, 'guarded residual unexpectedly subnormal'
        assert ref.exact(new[0]) + ref.exact(new[1]) == exact_product
        if expected is not None:
            assert new == expected[:2] and expected[2] == 1

start = time.monotonic()
before = [seal(p) for p in (Path(__file__).resolve(), REFERENCE, CORPUS, SOURCE)]
text = SOURCE.read_text()
original_body = ''' RPScalar ah,al,bh,bl;Split32(a,ah,al);Split32(b,bh,bl);hi=OrderedMultiply32(a,b);
 RPScalar e1=OrderedSubtract32(hi,OrderedMultiply32(ah,bh));
 RPScalar e2=OrderedSubtract32(e1,OrderedMultiply32(al,bh));
 RPScalar e3=OrderedSubtract32(e2,OrderedMultiply32(ah,bl));
 lo=OrderedSubtract32(OrderedMultiply32(al,bl),e3);'''
assert text.count(original_body) == 1
assert 'OrderedMultiply32(RPFromBits(0x45800800u),a)' in text
assert LOWER == ref.round_exact(ref.Fraction(1, 1 << 40))
assert UPPER == ref.round_exact(ref.Fraction(1 << 60))
raw = CORPUS.read_bytes()
magic, version, records = struct.unpack_from('<3I', raw)
assert (magic, version, records) == (ref.MAGIC, ref.VERSION, ref.RECORD_COUNT)
assert len(raw) == 12 + 24 * records
counts, families = Counter(), Counter()
residual_records = 0
for op, a, b, high, low, valid in struct.iter_unpack('<6I', raw[12:]):
    if op != 13:
        continue
    residual_records += 1
    check(a, b, 'existing-independent-residual-corpus', (high, low, valid))
assert residual_records == ref.RESIDUAL_COUNT == 53524

magnitudes = [LOWER-1, LOWER, LOWER+1, UPPER-1, UPPER, UPPER+1,
              0, 1, 0x00800000, 0x00800001, 0x3f000001,
              0x3f7fffff, 0x3f800000, 0x3f800001, 0x40000000,
              0x3f800400, 0x3ffff800]
for exponent in (-40, -39, -20, -1, 0, 1, 20, 40, 59):
    for fraction in (0, 1, 0x3ff, 0x400, 0x7ff, 0x800, 0x1000, 0x7fffff):
        magnitudes.append(((exponent+127) << 23) | fraction)
edges = sorted(set(magnitudes))
targeted = 0
for a in edges:
    for b in (LOWER-1, LOWER, LOWER+1, UPPER-1, UPPER, UPPER+1,
              0x3f7fffff, 0x3f800000, 0x3f800001, 0x40000000,
              0x3f800400, 0x3ffff800):
        for sign_a in (0, SIGN):
            for sign_b in (0, SIGN):
                check(a | sign_a, b | sign_b, 'guard-split-carry-tie-edges')
                targeted += 1

counterexamples = []
for a, b, category in [(0x00800001, 0x3f000001, 'underflow'),
                       (0x3f7fffff, 0x40000000, 'zero-residual-sign')]:
    old = original(a, b)
    fused = (old[0], fused_residual(a, b, old[0]))
    new, route = candidate(a, b)
    assert fused != old and new == old and route != 'guarded-nonzero'
    counterexamples.append({'a': f'{a:08x}', 'b': f'{b:08x}', 'category': category,
                            'original': [f'{x:08x}' for x in old],
                            'blanket_fma': [f'{x:08x}' for x in fused], 'route': route})
after = [seal(p) for p in (Path(__file__).resolve(), REFERENCE, CORPUS, SOURCE)]
assert before == after
assert all(counts[route] for route in ('guarded-nonzero', 'zero-residual-fallback',
                                      'outside-guard-fallback'))
result = {'pass': True, 'source_revision': 'b81e061fa081cdaed37eca3add950bd254bf296d',
          'device_workloads': 0, 'inputs': before, 'inputs_unchanged': True,
          'existing_residual_records': residual_records, 'targeted_records': targeted,
          'routes': dict(counts), 'families': dict(families), 'counterexamples': counterexamples,
          'wall_seconds': time.monotonic()-start,
          'contract': 'Guarded nonzero output words equal exact ordered original default; all other admitted inputs execute identical original fallback. Original five-min-subnormal enclosure retained.',
          'limitations': 'Finite Fraction model, not exhaustive or machine-checked equivalence; no actual compilation, shader selection/device execution, physical cost, rendering or qualification.'}
(WORK/'model-result.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k: result[k] for k in ('pass','existing_residual_records','targeted_records','routes','counterexamples','wall_seconds','limitations')},indent=2))
