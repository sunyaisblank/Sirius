"""Preserve the material value-interface compiler result without admission transfer."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import os
import subprocess
import zipfile

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-strict-ir-linking'
A = R / 'attestations/qualification/9cdafa03/native-cpu-strict-value-interface-compiler-gate'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,
                    sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not A.exists()
actions = json.loads((W / 'actions.json').read_text())
history = {}
for folder in sorted(W.iterdir()):
    path = folder / 'owner.json'
    if not folder.is_dir() or not path.exists():
        continue
    owner = json.loads(path.read_text())
    assert folder.name in actions and owner['status'] == 'terminal'
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    for birth in owner['observed_births']:
        try:
            current = (Path('/proc') / str(birth['pid']) / 'stat').read_text().rsplit(')', 1)[1].split()[19]
        except FileNotFoundError:
            continue
        assert current != birth['saved_start_ticks'], (folder.name, birth)
    subprocess.run(['git', 'cat-file', '-e', owner['source']['revision'] + '^{commit}'], check=True)
    history[folder.name] = dict(path=str(path.relative_to(R)), **pin(path), passed=owner['passed'],
        source=owner['source'], seconds=owner['wall_seconds'], peak_rss_kib=owner['peak_sampled_rss_kib'])
assert len(history) == 23
assert {k for k, v in history.items() if not v['passed']} == {'acquire', 'convert', 'value-lower', 'machine-dp'}
structural = json.loads((W / 'value-optimize-joins-result.json').read_text())
machine = json.loads((W / 'machine-features-result.json').read_text())
assert structural['adapter_calls'] == structural['packed_payload_allocas'] == structural['packed_tag_allocas'] == 0
assert structural['core_definitions'] == 8 and structural['core_calls'] == 1581
assert machine['no_FMA'] and len(machine['core_symbols']) == 8
for stage in ('transport', 'endpoint'):
    assert history['constrain-' + stage]['passed']
sources = [p for p in W.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
sources += [B / 'transport.ll', B / 'endpoint.ll']
sources += [R / ('bin/linux-gcc/native-cpu-inline-private-bank/ispc-source/' + s + '.ispc')
            for s in ('transport', 'endpoint', 'dopri_phase')]
sources += [R / 'out/native-cpu-inline-private-bank/packinghelper.isph',
            R / 'out/native-cpu-value-fact-stores/value_facts.isph',
            R / 'out/native-cpu-local-identities/identity_store.isph',
            R / 'out/native-cpu-exact-product-zero-residual/retained8.cpp']
products = [B / n for n in ('dp-value.bc', 'dp-value.joins.bc', 'dp-value.joins.O2.bc',
            'dp-value.features.bc', 'dp-value.features.o', 'transport.strict.bc', 'endpoint.strict.bc',
            'value_abi.so', 'value_abi.joins.so', 'constrain_remaining.so')]
A.mkdir(parents=True)
bundles = {}
for label, paths in (('source-and-receipts', sources), ('products', products)):
    entries = {str(p.relative_to(R)): pin(p) for p in sorted(set(paths))}
    target = A / (label + '.zip')
    with zipfile.ZipFile(target, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zipped:
        for name in entries:
            zipped.write(R / name, name)
    with zipfile.ZipFile(target) as zipped:
        assert set(zipped.namelist()) == set(entries)
        for name, expected in entries.items():
            with zipped.open(name) as stream:
                assert hashlib.file_digest(stream, 'sha256').hexdigest() == expected['sha256'], name
            assert zipped.getinfo(name).file_size == expected['bytes']
            assert pin(R / name) == expected, name
    assert target.stat().st_size < 100 * 1024 * 1024
    bundles[label] = dict(path=str(target.relative_to(R)), **pin(target), members=entries)
prior = R / 'attestations/qualification/97b3a761/native-cpu-strict-ir-linking-compiler-gate/MANIFEST.json'
manifest = dict(status='DP_PACKED_BANKS_ELIMINATED_NATIVE_OBJECT_EMITTED_ALL3_STRICT_IR_AVAILABLE',
    custody_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    source_authority='Each owner.source revision is immutable in this Git history. Current producers include explicit later fixes; no historical action or receipt is rewritten. Full prior DP input/link/O2 bytes are in the referenced earlier bundle.',
    previous_gate=dict(path=str(prior.relative_to(R)), **pin(prior)), owners=history,
    phase_passes=dict(Counter('passed' if v['passed'] else 'failed' for v in history.values())),
    DP=dict(value_output=structural['output'], optimized_IR=structural['ir'],
            adapter_calls=0, payload_banks=0, tag_banks=0, shared_cores=8, core_calls=1581,
            native=machine['output'], decoded_instructions=machine['decoded_instructions'],
            no_FMA=True, undefined_symbols=machine['undefined_symbols']),
    remaining_stages={s: json.loads((W / ('constrain-' + s + '-result.json')).read_text())
                      for s in ('transport', 'endpoint')},
    failure_dispositions=dict(acquire='Initial network timeout; exact immutable range continuation preserved.',
        convert='SDK opt omits SHA methods; only isolated SHA256 archive member linked.',
        value_lower='Identical local destination-pointer PHI normalized before bounds check; no differing-address PHI accepted.',
        machine_dp='Incorrect -fma target setting disabled dependent AVX512 features; removed only that feature token and llc setting, preserving strict arithmetic/no-contract and rejecting actual FMA.'),
    cleanup=json.loads((prior.parent / 'retirement.json').read_text()), bundles=bundles,
    decision='Material DP interface change reaches a real native object. Next integrate value interface for remaining two stages and qualify one complete coupled candidate; do not rerun the old failed frame or qualify partial stages as full acceptance.',
    limitations='Compiler/source/static native checks only. No new finite numerical/direct sticky-flag compatibility, full ISA/PC/FP-scope proof, C8/all32/ABBA/all8/provider/frame/hardware-domain/scientific/release admission. Existing public/GCC/raw/fallback objects and production source remain unchanged.')
write(A / 'MANIFEST.json', manifest)
print(json.dumps(dict(sealed=str(A.relative_to(R)), phases=len(history), phase_passes=manifest['phase_passes'],
                      bundles={k: v['bytes'] for k, v in bundles.items()}, DP=manifest['DP'])))
