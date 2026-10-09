"""Independent Fraction oracle for bounded high/residual product semantics.

Decode and nearest-grid rounding reuse the previously qualified Python rational
oracle; no C++ significand, lattice-subtraction or guard/sticky algorithm is used.
The new expectation is round(product), then round(exact product - decoded high).
Historical primitive corpus/expectations are read unchanged, never regenerated.
"""
from collections import Counter
import hashlib
import json
from pathlib import Path
import random
import runpy
import struct
import subprocess

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
HISTORY = ROOT.parent / "portable-arithmetic"
REFERENCE = runpy.run_path(str(HISTORY / "exact_reference.py"))
decode = REFERENCE["exact"]
round_exact = REFERENCE["round_exact"]
BOUND = 0x5d800000
SIGN = 0x80000000

def expected(a,b):
    if (a & 0x7fffffff) > BOUND or (b & 0x7fffffff) > BOUND:
        return 0,0,0
    product = decode(a) * decode(b)
    high = round_exact(product, (a ^ b) & SIGN)
    # RTE exact cancellation is +0, including signed-zero products minus high.
    low = round_exact(product - decode(high))
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

def corpus():
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

def main():
    for a,b,gold in ANALYTIC:
        assert expected(a,b) == gold, (hex(a),hex(b),expected(a,b),gold)
    historical_cases = [tuple(int(x,16) for x in line.split()) for line in (HISTORY/"corpus.txt").read_text().splitlines()]
    historical_gold = [(int(line,16),0,1) for line in (HISTORY/"expected.txt").read_text().splitlines()]
    assert len(historical_cases) == len(historical_gold) == 120744
    residual_cases = list(corpus())
    cases = historical_cases + [(13,a,b) for a,b,_ in residual_cases]
    gold = historical_gold + [expected(a,b) for a,b,_ in residual_cases]
    text = "".join(f"{op:x} {a:08x} {b:08x}\n" for op,a,b in cases)
    gold_text = "".join(f"{h:08x} {l:08x} {v:08x}\n" for h,l,v in gold)
    (ROOT/"corpus.txt").write_text(text)
    (ROOT/"expected.txt").write_text(gold_text)
    input_bytes = struct.pack("<I",len(cases)) + b"".join(struct.pack("<III",*c) for c in cases)
    expected_bytes = b"".join(struct.pack("<III",*g) for g in gold)
    (ROOT/"shader-input.bin").write_bytes(input_bytes)
    (ROOT/"shader-expected.bin").write_bytes(expected_bytes)
    process = subprocess.run([str(ROOT/"host_probe")],input=text,text=True,capture_output=True,check=True)
    actual = [tuple(int(x,16) for x in line.split()) for line in process.stdout.splitlines()]
    if len(actual) != len(gold):
        raise RuntimeError("complete host output count mismatch")
    failures = [{"row":i,"input":cases[i],"expected":gold[i],"actual":actual[i]} for i in range(len(gold)) if actual[i]!=gold[i]]
    lows = [g[1] for g in gold[len(historical_gold):] if g[2]]
    subnormal_lows = sum(0<(lo&0x7fffffff)<0x800000 for lo in lows)
    negative_zero_lows = lows.count(0x80000000)
    nonzero_lows = sum((lo&0x7fffffff)!=0 for lo in lows)
    assert subnormal_lows>0 and negative_zero_lows>0 and nonzero_lows>0
    report = {"scope":"canonical host promotion and bounded exact-product residual only; no GPU or production-route admission",
        "historical_primitive_cases":len(historical_cases),"new_residual_cases":len(residual_cases),"total_cases":len(cases),
        "analytic_selfchecks":len(ANALYTIC),"residual_categories":dict(Counter(c[2] for c in residual_cases)),
        "mismatches":len(failures),"failure_examples":failures[:20],
        "nonzero_rounded_residuals":nonzero_lows,"negative_zero_residuals":negative_zero_lows,
        "FTZ_low_mutant_mismatches":subnormal_lows,"force_positive_zero_low_mutant_mismatches":sum(lo!=0 for lo in lows),
        "header_sha256":hashlib.sha256((REPO/"src/sirius/kernels/portable_binary32.h").read_bytes()).hexdigest(),
        "oracle_source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "reused_rational_oracle_sha256":hashlib.sha256((HISTORY/"exact_reference.py").read_bytes()).hexdigest(),
        "input_sha256":hashlib.sha256(input_bytes).hexdigest(),"expected_sha256":hashlib.sha256(expected_bytes).hexdigest()}
    (ROOT/"host-report.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    raise SystemExit(bool(failures))

if __name__ == "__main__":
    main()
