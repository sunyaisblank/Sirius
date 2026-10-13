"""Generate independent exact binary32 expectations, never native float arithmetic.

Binary words decode to Fraction. Quotient/remainder rounds to the nearest binary32
lattice with ties to even; sqrt independently compares exact squared midpoints.
For bounded products, low is round(exact product - decode(round(exact product))).
The significand/guard/sticky and two-word subtraction implementation is not used.

Fixture v1: little-endian uint32 magic("PB32"), version, record count, followed by
six uint32 words per record: op, a, b, expected high, expected low, expected valid.
Ops 0..12 cover scalar arithmetic/comparisons/conversions/bit operations; op13
covers finite product residuals for |a|,|b| <= 2^60. Arithmetic NaNs use canonical
0x7fc00000, bit operations preserve payloads, and exception flags are not modeled.
This finite corpus is not a proof for every binary32 operand combination.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import random
import struct
import tempfile

MAGIC = 0x32334250
VERSION = 1
PRIMITIVE_COUNT = 120744
RESIDUAL_COUNT = 53524
RECORD_COUNT = PRIMITIVE_COUNT + RESIDUAL_COUNT

SIGN = 0x80000000
MASK = 0x7fffffff
INF = 0x7f800000
NAN = 0x7fc00000

def classify(a):
    a &= MASK
    return 0 if a == 0 else 1 if a < 0x800000 else 2 if a < INF else 3 if a == INF else 4

def exact(a):
    m = a & MASK
    e = m >> 23
    f = m & 0x7fffff
    if e == 255:
        return None
    sig = f if e == 0 else f | 0x800000
    exponent = -149 if e == 0 else e - 150
    value = Fraction(sig << exponent, 1) if exponent >= 0 else Fraction(sig, 1 << -exponent)
    return -value if a & SIGN else value

def exponent_of(x):
    e = x.numerator.bit_length() - x.denominator.bit_length()
    if e >= 0:
        if x.numerator < x.denominator << e:
            e -= 1
    elif x.numerator << -e < x.denominator:
        e -= 1
    return e

def pack_rounded(sign, e, q):
    if q == 0:
        return sign
    if e < -126:
        if q < 0x800000:
            return sign | q
        e = -126
    if q >= 0x1000000:
        q >>= 1
        e += 1
    if e > 127:
        return sign | INF
    return sign | ((e + 127) << 23) | (q & 0x7fffff)

def round_exact(x, zero_sign=0):
    if x == 0:
        return zero_sign
    sign = SIGN if x < 0 else 0
    x = abs(x)
    e = exponent_of(x)
    unit_exponent = max(e - 23, -149)
    n, d = x.numerator, x.denominator
    if unit_exponent < 0:
        n <<= -unit_exponent
    else:
        d <<= unit_exponent
    q, r = divmod(n, d)
    if 2 * r > d or (2 * r == d and q & 1):
        q += 1
    return pack_rounded(sign, e, q)

def round_sqrt(x):
    e = exponent_of(x) // 2
    unit_exponent = max(e - 23, -149)
    n, d = x.numerator, x.denominator
    if unit_exponent < 0:
        n <<= -2 * unit_exponent
    else:
        d <<= 2 * unit_exponent
    q = math.isqrt(n // d)
    midpoint_test = 4 * n - d * (2 * q + 1) ** 2
    if midpoint_test > 0 or (midpoint_test == 0 and q & 1):
        q += 1
    return pack_rounded(0, e, q)

def add(a, b):
    ca, cb = classify(a), classify(b)
    if ca == 4 or cb == 4:
        return NAN
    if ca == 3 or cb == 3:
        return NAN if ca == cb == 3 and (a ^ b) & SIGN else a if ca == 3 else b
    return round_exact(exact(a) + exact(b), (a & b & SIGN) if ca == cb == 0 else 0)

def compare(a, b):
    if classify(a) == 4 or classify(b) == 4:
        return None
    if classify(a) == 3:
        if classify(b) == 3:
            return 0 if a == b else -1 if a & SIGN else 1
        return -1 if a & SIGN else 1
    if classify(b) == 3:
        return 1 if b & SIGN else -1
    return (exact(a) > exact(b)) - (exact(a) < exact(b))

def scalar_expected(op, a, b):
    ca, cb = classify(a), classify(b)
    if op == 0:
        return add(a, b)
    if op == 1:
        return add(a, b ^ SIGN)
    if op == 2:
        sign = (a ^ b) & SIGN
        if ca == 4 or cb == 4 or (ca == 3 and cb == 0) or (cb == 3 and ca == 0):
            return NAN
        if ca == 3 or cb == 3:
            return sign | INF
        return round_exact(exact(a) * exact(b), sign)
    if op == 3:
        sign = (a ^ b) & SIGN
        if ca == 4 or cb == 4 or ca == cb == 3 or ca == cb == 0:
            return NAN
        if ca == 3 or cb == 0:
            return sign | INF
        if ca == 0 or cb == 3:
            return sign
        return round_exact(exact(a) / exact(b), sign)
    if op == 4:
        if ca == 4 or a & SIGN and ca != 0:
            return NAN
        return a if ca in (0, 3) else round_sqrt(exact(a))
    if op in (5, 6, 7):
        c = compare(a, b)
        return int(c is not None and (c == 0 if op == 5 else c < 0 if op == 6 else c <= 0))
    if op == 8:
        return round_exact(Fraction(a))
    if op == 9:
        return round_exact(Fraction(a - (1 << 32) if a & SIGN else a))
    if op == 10:
        return a & MASK
    if op == 11:
        return a ^ SIGN
    if op == 12:
        return ca
    raise ValueError(op)

def primitive_cases():
    edge_magnitudes = [0, 1, 2, 3, 4, 7, 8, 0x3fffff, 0x7ffffe, 0x7fffff,
                       0x800000, 0x800001, 0x800002, 0x3eaaaaab, 0x3f000000,
                       0x3f7ffffe, 0x3f7fffff, 0x3f800000, 0x3f800001,
                       0x3f800002, 0x3fffffff, 0x40000000, 0x40400000,
                       0x40800000, 0x41100000, 0x4b000000, 0x4b7fffff,
                       0x4b800000, 0x4b800001, 0x5d800000, 0x7b800000,
                       0x7f7ffffe, 0x7f7fffff, INF, 0x7f800001, NAN, 0x7fffffff]
    edges = [x | sign for x in edge_magnitudes for sign in (0, SIGN)]
    rng = random.Random(0x5349524955530048)
    for op in (0, 1, 2, 3, 5, 6, 7):
        for a in edges:
            for b in edges:
                yield op, a, b, "edge-cross"
        for _ in range(4096):
            yield op, rng.getrandbits(32), rng.getrandbits(32), "raw-fixed-seed"
        for _ in range(4096):
            a = rng.randrange(1, 0x800000) | rng.choice((0, SIGN))
            b = rng.randrange(1, 0x800000) | rng.choice((0, SIGN))
            yield op, a, b, "subnormal-fixed-seed"
    for op in (4, 8, 9, 10, 11, 12):
        for a in edges + [0xffffffff, 0x7fffffff, 0x01000001, 0x01000003]:
            yield op, a, 0, "unary-edge"
        for _ in range(4096):
            yield op, rng.getrandbits(32), 0, "unary-fixed-seed"
    # Exact halfway witnesses and transitions; RTE and FTZ mutations differ here.
    for a, b in [(0x3f800000, 0x33800000), (0x3f800001, 0x33800000),
                 (1, 0x3f000000), (3, 0x3f000000),
                 (0x7fffff, 1), (0x800000, 0x7fffff)]:
        for op in range(4):
            yield op, a, b, "declared-halfway-boundary"

BOUND = 0x5d800000

def product_expected(a,b):
    if (a & 0x7fffffff) > BOUND or (b & 0x7fffffff) > BOUND:
        return 0,0,0
    product = exact(a) * exact(b)
    high = round_exact(product, (a ^ b) & SIGN)
    # RTE exact cancellation is +0, including signed-zero products minus high.
    low = round_exact(product - exact(high))
    return high,low,1

ANALYTIC = [
    (0x3f800001,0x3f800001,(0x3f800002,0x28800000,1)),
    (0xbf800001,0x3f800001,(0xbf800002,0xa8800000,1)),
    (0x3f800001,0x3fc00000,(0x3fc00002,0xb3800000,1)),
    (0xbf800001,0x3fc00000,(0xbfc00002,0x33800000,1)),
    (0xbf800000,0x3f800000,(0xbf800000,0,1)),
    (0x80000000,0x3f800000,(0x80000000,0,1)),
    (0x00000001,0x3f000000,(0,0,1)),
    (0x80000001,0x3f000000,(0x80000000,0x80000000,1)),
    (0x00800000,0x34000000,(1,0,1)),
    (0x007fffff,0x34000000,(1,0x80000000,1)),
    (0x00800001,0x33800000,(1,0x80000000,1)),
    (0x80800001,0x33800000,(0x80000001,0,1)),
    (0x17800001,0x35800001,(0x0d800002,8,1)),
    (0x3f800400,0x3ffff800,(0x40000000,0xb3000000,1)),
    (BOUND,BOUND,(0x7b800000,0,1)),
    (BOUND+1,0x3f800000,(0,0,0)),
]

def residual_cases():
    magnitudes = [0,1,2,3,4,0x3fffff,0x7ffffe,0x7fffff,0x800000,0x800001,
        0x800002,0x17800000,0x17800001,0x28800000,0x33800000,0x34000000,
        0x35800001,0x3f000000,0x3f7fffff,0x3f800000,0x3f800001,0x3f800002,
        0x3f800400,0x3fc00000,0x3ffff800,0x3fffffff,0x40000000,
        0x4b800000,BOUND-1,BOUND,BOUND+1,0x7f800000,0x7fc00000]
    edges = [a|s for a in magnitudes for s in (0,SIGN)]
    for a in edges:
        for b in edges:
            yield a,b,"edge-cross"
    rng = random.Random(0x5349524955530049)
    def signed(word):
        return word | rng.choice((0,SIGN))
    for _ in range(16384):
        yield signed(rng.randrange(BOUND+1)),signed(rng.randrange(BOUND+1)),"bounded-raw-fixed-seed"
    for _ in range(8192):
        yield signed(rng.randrange(1,0x800000)),signed(rng.randrange(BOUND+1)),"subnormal-input-fixed-seed"
    exponents = (-126,-125,-100,-80,-60,-40,-20,-1,0,1,20,40,59)
    for _ in range(16384):
        a = ((rng.choice(exponents)+127)<<23) | rng.randrange(0x800000)
        b = ((rng.choice(exponents)+127)<<23) | rng.randrange(0x800000)
        yield signed(a),signed(b),"scaled-mantissa-fixed-seed"
    for _ in range(8192):
        total = rng.choice((-152,-151,-150,-149,-148,-127,-126,-125,-103,-100,-80,0,60,119))
        if total == 119:
            # At exponent 60 the bounded operand must be exactly 2^60.
            yield signed(BOUND),signed((186<<23)|rng.randrange(0x800000)),"product-rounding-boundaries"
            continue
        ea = rng.randrange(max(-126,total-59),min(59,total+126)+1)
        eb = total-ea
        a = ((ea+127)<<23) | rng.randrange(0x800000)
        b = ((eb+127)<<23) | rng.randrange(0x800000)
        yield signed(a),signed(b),"product-rounding-boundaries"
    for a,b,_ in ANALYTIC:
        yield a,b,"analytic-selfcheck"

def selfcheck():
    scalar = [(0,0x3f800000,0x33800000,0x3f800000),
              (0,0x3f800001,0x33800000,0x3f800002),
              (2,1,0x3f000000,0),(2,3,0x3f000000,2),
              (2,0x80000001,0x3f000000,0x80000000),
              (0,0x007fffff,1,0x00800000),(1,0x00800000,0x007fffff,1),
              (8,0x01000003,0,0x4b800002),(4,0x40800000,0,0x40000000),
              (4,0x80000000,0,0x80000000)]
    for op,a,b,gold in scalar:
        if scalar_expected(op,a,b) != gold:
            raise RuntimeError("analytic scalar oracle selfcheck failed")
    for a,b,gold in ANALYTIC:
        if product_expected(a,b) != gold:
            raise RuntimeError("analytic residual oracle selfcheck failed")

def records():
    for op,a,b,category in primitive_cases():
        yield op,a,b,scalar_expected(op,a,b),0,1,"primitive/"+category
    for a,b,category in residual_cases():
        high,low,valid = product_expected(a,b)
        yield 13,a,b,high,low,valid,"residual/"+category

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True,type=Path,help="generated binary fixture path")
    args = parser.parse_args()
    selfcheck()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    header = struct.pack("<III",MAGIC,VERSION,RECORD_COUNT)
    digest = hashlib.sha256(header)
    categories = Counter()
    primitive_count = residual_count = negative_zero_lows = subnormal_lows = nonzero_lows = 0
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(prefix=args.output.name+".",dir=args.output.parent,delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(header)
            for row,record in enumerate(records()):
                if row >= RECORD_COUNT:
                    raise RuntimeError("independent corpus exceeds versioned finite record count")
                words = record[:6]
                packed = struct.pack("<IIIIII",*words)
                stream.write(packed)
                digest.update(packed)
                categories[record[6]] += 1
                if words[0] != 13:
                    primitive_count += 1
                else:
                    residual_count += 1
                    if words[5]:
                        low = words[4]
                        negative_zero_lows += low == SIGN
                        subnormal_lows += 0 < (low&MASK) < 0x800000
                        nonzero_lows += (low&MASK) != 0
            if (primitive_count,residual_count) != (PRIMITIVE_COUNT,RESIDUAL_COUNT):
                raise RuntimeError("independent corpus is incomplete")
            if not negative_zero_lows or not subnormal_lows or not nonzero_lows:
                raise RuntimeError("residual hard regimes are empty")
            stream.flush()
        temporary.replace(args.output)
    except BaseException:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise
    print(json.dumps({"magic":"PB32","version":VERSION,"records":RECORD_COUNT,
        "primitive_cases":primitive_count,"residual_cases":residual_count,
        "negative_zero_residuals":negative_zero_lows,"subnormal_residuals":subnormal_lows,
        "nonzero_rounded_residuals":nonzero_lows,"categories":dict(categories),
        "bytes":args.output.stat().st_size,"fixture_sha256":digest.hexdigest()},sort_keys=True))

if __name__ == "__main__":
    main()
