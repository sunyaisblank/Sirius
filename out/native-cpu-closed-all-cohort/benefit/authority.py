"""Pin the exercised original AU benefit chain; change current B binding/scope only."""
from pathlib import Path
import hashlib
_AU_BENEFIT_SOURCE=Path.cwd()/'out/native-cpu-shared-value-cohort/benefit'
_AU_BENEFIT_PINS={'authority.py': '898ec64ce1e9cc68f43a0bb158332cb3493fbd96e1e5ac585d7437787abb5c71', 'correctness.py': '4a4c0a64767accdfd2a68bd0e51e409925b0a5dd611a7f65cf72c75f73b5effc', 'prepare.py': '7b46dc99e762626f40df8d22991b16a30f852c8167f9fb98929211b770a63c3b', 'readers/decimal_check.py': '37705aa878ac58daa25d6f14588bdc06be0dbec8a825e0aa59103844aa7f3928', 'readers/native_evidence.py': '76dbffcce36e1d4a11161085631e06be1aa0e561046ca4b404201b624a22c453', 'readers/numeric_gates.py': '78c899fde5b7846e28c857323fa68af8fc70e8c6a24097e404a70cbb9acbdb2b', 'readers/read_results.py': '57df180c7f2dfec81e0f8e2f157cff5657d12149779ba840d1d4e24056611c18', 'readers/run_abba.py': '5dfa4d342239391c54a88eb592d3e97b0f9ad3216eeb6fc47b7d82764ee365ca', 'runtime/owned_guard.py': '19f2ba5d29da31c8be60d3644730c0b634288565f212fcafa872354d0004aaee', 'runtime/run_case.py': 'e451025d1518eb558adce27d9369d1f27ebd4aa9bc03af818856cfb98e5d1086', 'runtime/run_owned.py': '22503c90536d4f93fb1f8c22a3a6c39a9b8237f7a1099013346bafd9a24a907f', 'readers/expected-readback-layout.json': 'a40ce35b744ff05c987833de8337ec6f5445daea9553348cdb9209ad78314392'}
for _name,_digest in _AU_BENEFIT_PINS.items():
    assert hashlib.sha256((_AU_BENEFIT_SOURCE/_name).read_bytes()).hexdigest()==_digest,_name
assert hashlib.sha256((Path(__file__).resolve().parent/'readers/expected-readback-layout.json').read_bytes()).hexdigest()==_AU_BENEFIT_PINS['readers/expected-readback-layout.json']
exec(compile((_AU_BENEFIT_SOURCE/'authority.py').read_text(),str(_AU_BENEFIT_SOURCE/'authority.py'),'exec'))
_previous_AX_load=load

def load(name):
    text=_previous_AX_load(name)
    return text.replace('out/native-cpu-shared-value-cohort','out/native-cpu-closed-all-cohort').replace('out/native-cpu-shared-value-provider/provider','out/native-cpu-closed-all-provider/provider')

_previous_AX_authorities=authority_files
def authority_files():
    return _previous_AX_authorities()+[_AU_BENEFIT_SOURCE/name for name in _AU_BENEFIT_PINS]

_previous_AX_prepare_source=prepare_source
def prepare_source(text):
    text=_previous_AX_prepare_source(text)
    return text.replace('current H host/AU','current H host/AX')

_previous_AX_reader_scope=reader_scope
def reader_scope(text):
    text=_previous_AX_reader_scope(text)
    before='Current AU preserves AT public/scalar/admission and AS public raw/math fallbacks; only positive-job v4 private workers use the new strict raw3/shared-nine-core value interface inside the full saved/masked/restored FP scope. Current genuine-job C8/all32 is its first numerical gate on the same unchanged H host. Original seven controls and final all8 receipt remain required after fresh22+Decimal. Cold-frame/scientific/release gates remain separate.'
    after='Current AX preserves AT public/scalar/admission and AS public raw/GCC/v3/AVX2 fallbacks; genuine positive-job v4 private workers select closed canonical T/E/DP with one existing strict local-effects regcall math9 inside the complete saved/masked/restored FP scope. The current original C8/all32 passed on the unchanged H host: Transport12 selected original scalar work, while Endpoint and DopriPhase supplied actual closed-raw/local-core PC coverage. No Transport SIMD execution is inferred from that C8. This sole original ABBA must satisfy all22 numeric/full-buffer gates and independent Decimal agreement before original seven controls/final all8/cold-frame/scientific/release gates; no earlier benefit or frame acceptance transfers.'
    assert text.count(before)==1
    return text.replace(before,after,1)
