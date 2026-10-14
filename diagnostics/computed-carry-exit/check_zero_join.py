"""Finite independent Fraction-oracle check of the skipped zero-addend join."""
from pathlib import Path
import hashlib
import importlib.util
import itertools
import json
import sys

sys.dont_write_bytecode = True
root = Path.cwd().resolve()
work = Path(__file__).resolve().parent
reference = root / 'tests/support/portable_binary32/reference.py'
spec = importlib.util.spec_from_file_location('issue101_fraction_reference', reference)
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
values = set()
for exponent, mantissa, sign in itertools.product(range(255),
        (0, 1, 0x3fffff, 0x7ffffe, 0x7fffff), (0, 0x80000000)):
    values.add(sign | (exponent << 23) | mantissa)
count = 0
for limb, carry in itertools.product(sorted(values), (0, 0x80000000)):
    a, b = limb, carry
    if (a & 0x7fffffff) < (b & 0x7fffffff):
        a, b = b, a
    high = ref.add(a, b)
    displacement = ref.add(high, a ^ 0x80000000)
    low = ref.add(b, displacement ^ 0x80000000)
    expected = (limb if limb & 0x7fffffff else limb & carry, carry)
    assert (high, low) == expected, (limb, carry, high, low, expected)
    assert displacement == 0, (limb, carry, displacement)
    count += 1
special = []
for limb, carry in itertools.product(
        (0x7f800000, 0xff800000, 0x7f800001, 0x7fc00000, 0xffffffff), (0, 0x80000000)):
    high = ref.add(limb, carry)
    displacement = ref.add(high, limb ^ 0x80000000)
    low = ref.add(carry, displacement ^ 0x80000000)
    assert (high & 0x7fffffff) > 0x7d800000
    special.append([hex(limb), hex(carry), hex(high), hex(low)])
result = {'finite_word_pairs': count,
    'plain_three_add_chain_matches_normal_zero_arm': True,
    'displacement_is_positive_zero': True,
    'special_unvisited_pairs_refuse_private_bound': len(special), 'special_words': special,
    'reference_sha256': hashlib.sha256(reference.read_bytes()).hexdigest(),
    'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'scope': 'Independent Fraction-oracle zero-join checks only; no candidate shader execution, all-operands proof, compiler or runtime verdict.'}
(work / 'zero-join-check.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
